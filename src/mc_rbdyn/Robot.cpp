/*
 * Copyright 2015-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_rbdyn/Robot.h>
#include <mc_rbdyn/RobotModule.h>
#include <mc_rbdyn/Robots.h>
#include <mc_rbdyn/SCHAddon.h>
#include <mc_rbdyn/Surface.h>
#include <mc_rbdyn/ZMP.h>
#include <mc_rbdyn/surface_utils.h>
#include <mc_rtc/constants.h>
#include <mc_rtc/logging.h>
#include <mc_rtc/pragma.h>

#include <RBDyn/CoM.h>
#include <RBDyn/EulerIntegration.h>
#include <RBDyn/FA.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>

#include <sch/S_Object/S_Cylinder.h>
#include <sch/S_Object/S_Superellipsoid.h>

#include <boost/filesystem.hpp>
namespace bfs = boost::filesystem;

#include <fstream>
#include <tuple>

namespace
{

using bound_t = std::vector<std::vector<double>>;
using bounds_t = std::tuple<bound_t, bound_t, bound_t, bound_t, bound_t, bound_t>;
using accelerationBounds_t = std::tuple<bound_t, bound_t>;
using jerkBounds_t = std::tuple<bound_t, bound_t>;
using torqueDerivativeBounds_t = std::tuple<bound_t, bound_t>;
using rm_bounds_t = mc_rbdyn::RobotModule::bounds_t;
using rm_bound_t = rm_bounds_t::value_type;

using jt_method = int (rbd::Joint::*)() const;

bound_t fill_bound(const rbd::MultiBody & mb,
                   const std::string & name,
                   const rm_bound_t & bound_in,
                   jt_method def_size,
                   double def_value,
                   double ff_def_value)
{
  bound_t res;
  res.reserve(static_cast<size_t>(mb.nrJoints()));
  for(const auto & j : mb.joints())
  {
    res.emplace_back(((j).*(def_size))(), j.type() == rbd::Joint::Free ? ff_def_value : def_value);
    if(bound_in.count(j.name()))
    {
      const auto & b_ref = bound_in.at(j.name());
      auto & b = res.back();
      if(b_ref.size() != b.size())
      {
        mc_rtc::log::error_and_throw("{} provided bound size ({}) different from expected size ({}) for joint {}", name,
                                     b_ref.size(), b.size(), j.name());
      }
      res.back() = bound_in.at(j.name());
    }
  }
  return res;
}

/** Generate bounds compatible with the given MultiBody
 *
 * If bounds is provided, use the values provided to build the bounds.
 *
 * Otherwise, default bounds are returned.
 */
bounds_t bounds(const rbd::MultiBody & mb, const rm_bounds_t & bounds)
{
  return std::make_tuple(fill_bound(mb, "lower position", bounds.at(0), &rbd::Joint::params, -INFINITY, -INFINITY),
                         fill_bound(mb, "upper position", bounds.at(1), &rbd::Joint::params, INFINITY, INFINITY),
                         fill_bound(mb, "lower velocity", bounds.at(2), &rbd::Joint::dof, -INFINITY, -INFINITY),
                         fill_bound(mb, "upper velocity", bounds.at(3), &rbd::Joint::dof, INFINITY, INFINITY),
                         fill_bound(mb, "lower torque", bounds.at(4), &rbd::Joint::dof, -INFINITY, 0),
                         fill_bound(mb, "upper torque", bounds.at(5), &rbd::Joint::dof, INFINITY, 0));
}

/** Generate acceleration bounds compatible with the given MultiBody
 *
 * If bounds is provided, use the values provided to build the bounds.
 *
 * Otherwise, default bounds are returned.
 */
accelerationBounds_t acceleration_bounds(const rbd::MultiBody & mb, const rm_bounds_t & bounds)
{
  rm_bound_t default_bound = {};
  auto safe_bounds = [&bounds, &default_bound](size_t idx) -> const rm_bound_t & {
    if(idx < bounds.size())
    {
      return bounds[idx];
    }
    return default_bound;
  };
  return std::make_tuple(fill_bound(mb, "lower acceleration", safe_bounds(0), &rbd::Joint::dof, -INFINITY, -INFINITY),
                         fill_bound(mb, "upper acceleration", safe_bounds(1), &rbd::Joint::dof, INFINITY, INFINITY));
}

/** Generate jerk bounds compatible with the given MultiBody
 *
 * If bounds is provided, use the values provided to build the bounds.
 *
 * Otherwise, default bounds are returned.
 */
jerkBounds_t jerk_bounds(const rbd::MultiBody & mb, const rm_bounds_t & bounds)
{
  rm_bound_t default_bound = {};
  auto safe_bounds = [&bounds, &default_bound](size_t idx) -> const rm_bound_t & {
    if(idx < bounds.size())
    {
      return bounds[idx];
    }
    return default_bound;
  };
  return std::make_tuple(fill_bound(mb, "lower jerk", safe_bounds(0), &rbd::Joint::dof, -INFINITY, -INFINITY),
                         fill_bound(mb, "upper jerk", safe_bounds(1), &rbd::Joint::dof, INFINITY, INFINITY));
}

/** Generate torque-derivative bounds compatible with the given MultiBody
 *
 * If bounds is provided, use the values provided to build the bounds.
 *
 * Otherwise, default bounds are returned.
 */
torqueDerivativeBounds_t torqueDerivative_bounds(const rbd::MultiBody & mb, const rm_bounds_t & bounds)
{
  rm_bound_t default_bound = {};
  auto safe_bounds = [&bounds, &default_bound](size_t idx) -> const rm_bound_t & {
    if(idx < bounds.size())
    {
      return bounds[idx];
    }
    return default_bound;
  };
  return std::make_tuple(
      fill_bound(mb, "lower torque-derivative", safe_bounds(0), &rbd::Joint::dof, -INFINITY, -INFINITY),
      fill_bound(mb, "upper torque-derivative", safe_bounds(1), &rbd::Joint::dof, INFINITY, INFINITY));
}

template<typename schT, typename mapT>
void loadSCH(const mc_rbdyn::Robot & robot,
             const std::map<std::string, std::pair<std::string, std::string>> & urls,
             schT * (*sch_load_fn)(const std::string &),
             mapT & data_,
             std::map<std::string, sva::PTransformd> & cTfs)
{
  for(const auto & cH : urls)
  {
    const std::string & cHName = cH.first;
    const std::string & parent = cH.second.first;
    const std::string & cHURL = cH.second.second;
    if(robot.hasBody(parent))
    {
      auto poly = std::shared_ptr<schT>(sch_load_fn(cHURL));
      sch::mc_rbdyn::transform(*poly, robot.bodyPosW()[robot.bodyIndexByName(parent)]);
      data_[cHName] = {parent, poly};
      cTfs[cHName] = sva::PTransformd::Identity();
    }
  }
}

template<typename mapT>
void fixSCH(const mc_rbdyn::Robot & robot, mapT & data_, const std::map<std::string, sva::PTransformd> & tfs)
{
  for(const auto & d : data_)
  {
    const auto & pos = robot.bodyPosW(d.second.first);
    if(tfs.count(d.first))
    {
      sch::mc_rbdyn::transform(*d.second.second, tfs.at(d.first) * pos);
    }
    else
    {
      sch::mc_rbdyn::transform(*d.second.second, pos);
    }
  }
}

bool VisualToConvex(const std::string & robot,
                    const std::string & cName,
                    const std::string & bName,
                    const rbd::parsers::Visual & visual,
                    std::map<std::string, mc_rbdyn::Robot::convex_pair_t> & convexes,
                    std::map<std::string, sva::PTransformd> & collisionTransforms)
{
  // Ignore visual types that we cannot easily map to SCH
  if(visual.geometry.type == rbd::parsers::Geometry::Type::UNKNOWN
     || visual.geometry.type == rbd::parsers::Geometry::Type::MESH)
  {
    return false;
  }
  // If we already have a convex with the same name, discard loading
  if(convexes.count(cName) != 0)
  {
    mc_rtc::log::warning("While loading {}, a convex was already provided for collision geometry specified in URDF",
                         robot);
    return false;
  }
  auto fromBox = [&]() {
    const auto & box = boost::get<rbd::parsers::Geometry::Box>(visual.geometry.data);
    convexes[cName] = {bName, std::make_shared<sch::S_Box>(box.size.x(), box.size.y(), box.size.z())};
  };
  auto fromCylinder = [&]() {
    const auto & cyl = boost::get<rbd::parsers::Geometry::Cylinder>(visual.geometry.data);
    convexes[cName] = {bName, std::make_shared<sch::S_Cylinder>(sch::Point3(0, 0, -cyl.length / 2),
                                                                sch::Point3(0, 0, cyl.length / 2), cyl.radius)};
  };
  auto fromSphere = [&]() {
    const auto & sph = boost::get<rbd::parsers::Geometry::Sphere>(visual.geometry.data);
    convexes[cName] = {bName, std::make_shared<sch::S_Sphere>(sph.radius)};
  };
  auto fromSuperEllipsoid = [&]() {
    const auto & sel = boost::get<rbd::parsers::Geometry::Superellipsoid>(visual.geometry.data);
    convexes[cName] = {bName, std::make_shared<sch::S_Superellipsoid>(sel.size.x(), sel.size.y(), sel.size.z(),
                                                                      sel.epsilon1, sel.epsilon2)};
  };
  switch(visual.geometry.type)
  {
    case rbd::parsers::Geometry::Type::BOX:
      fromBox();
      break;
    case rbd::parsers::Geometry::Type::CYLINDER:
      fromCylinder();
      break;
    case rbd::parsers::Geometry::Type::SPHERE:
      fromSphere();
      break;
    case rbd::parsers::Geometry::Type::SUPERELLIPSOID:
      fromSuperEllipsoid();
      break;
    default:
      return false;
  }
  collisionTransforms[cName] = visual.origin;
  return true;
}

} // namespace

namespace mc_rbdyn
{

// We can safely ignore those since they are due to different index types and
// our index never go near unsafe territories
MC_RTC_diagnostic_push
MC_RTC_diagnostic_ignored(GCC, "-Wsign-conversion", ClangOnly, "-Wshorten-64-to-32")

Robot::Robot(const std::string & name,
             Robots & robots,
             unsigned int robots_idx,
             bool loadFiles,
             const sva::PTransformd * base,
             const std::string & bName)
: robots_(&robots), robots_idx_(robots_idx), name_(name)
{
  const auto & module_ = module();

  if(base)
  {
    std::string baseName = bName.empty() ? mb().body(0).name() : bName;
    mb() = mbg().makeMultiBody(baseName, mb().joint(0).type() == rbd::Joint::Fixed, *base);
    mbc() = rbd::MultiBodyConfig(mb());
  }

  mbc().gravity = mc_rtc::constants::gravity;
  mbc().zero(mb());
  {
    auto initQ = mbc().q;
    const auto & stance = module_.stance();
    for(size_t i = 0; i < mb().joints().size(); ++i)
    {
      const auto & j = mb().joint(static_cast<int>(i));
      if(stance.count(j.name()))
      {
        const auto & jQ = stance.at(j.name());
        if(initQ[i].size() != jQ.size())
        {
          mc_rtc::log::error_and_throw(
              "Missmatch between RobotModule stance for joint {}\nStance provides {} values but should be {}", j.name(),
              jQ.size(), initQ[i].size());
        }
        initQ[i] = jQ;
      }
    }
    if(initQ[0].size())
    {
      const auto & attitude = module_.default_attitude();
      initQ[0] = {std::begin(attitude), std::end(attitude)};
    }
    mbc().q = initQ;
    forwardKinematics();
  }

  bodyTransforms_.resize(mb().bodies().size());
  const auto & bbts =
      base ? mbg().bodiesBaseTransform(mb().body(0).name(), *base) : mbg().bodiesBaseTransform(mb().body(0).name());
  for(size_t i = 0; i < mb().bodies().size(); ++i)
  {
    const auto & b = mb().body(static_cast<int>(i));
    bodyTransforms_[i] = bbts.at(b.name());
  }

  if(module_.bounds().size() != 6)
  {
    mc_rtc::log::error_and_throw<std::invalid_argument>("The (urdf)-bounds of RobotModule \"{}\" have a size of {} "
                                                        "instead of 6 (ql, qu, vl, vu, tl, tu).",
                                                        module_.name, module_.bounds().size());
  }
  std::tie(ql_, qu_, vl_, vu_, tl_, tu_) = bounds(mb(), module_.bounds());

  if(module_.accelerationBounds().size() != 0 && module_.accelerationBounds().size() != 2)
  {
    mc_rtc::log::error_and_throw<std::invalid_argument>(
        "The additional acceleration bounds of RobotModule \"{}\" have a size of {} "
        "instead of 2 ([al, au]).",
        module_.name, module_.accelerationBounds().size());
  }
  std::tie(al_, au_) = acceleration_bounds(mb(), module_.accelerationBounds());

  if(module_.jerkBounds().size() != 0 && module_.jerkBounds().size() != 2)
  {
    mc_rtc::log::error_and_throw<std::invalid_argument>(
        "The additional jerk bounds of RobotModule \"{}\" have a size of {} "
        "instead of 2 ([jl, ju]).",
        module_.name, module_.jerkBounds().size());
  }
  std::tie(jl_, ju_) = jerk_bounds(mb(), module_.jerkBounds());

  if(module_.torqueDerivativeBounds().size() != 0 && module_.torqueDerivativeBounds().size() != 2)
  {
    mc_rtc::log::error_and_throw<std::invalid_argument>(
        "The additional torque-derivative bounds of RobotModule \"{}\" have a size of {} "
        "instead of 2 ([tdl, tdu]).",
        module_.name, module_.torqueDerivativeBounds().size());
  }
  std::tie(tdl_, tdu_) = torqueDerivative_bounds(mb(), module_.torqueDerivativeBounds());

  if(loadFiles)
  {
    loadSCH(*this, module_.convexHull(), &sch::mc_rbdyn::Polyhedron, convexes_, collisionTransforms_);
    for(const auto & c : module_._collision)
    {
      const auto & body = c.first;
      const auto & collisions = c.second;
      if(collisions.size() == 1)
      {
        VisualToConvex(name_, body, body, collisions[0], convexes_, collisionTransforms_);
        continue;
      }
      size_t added = 0;
      for(const auto & col : collisions)
      {
        if(VisualToConvex(name_, body + "_" + std::to_string(added), body, col, convexes_, collisionTransforms_))
        {
          added++;
        }
      }
    }
    for(const auto & o : module_.collisionObjects())
    {
      if(convexes_.count(o.first) != 0)
      {
        mc_rtc::log::warning("While loading {}, another object named {} was already loaded, the object specified in "
                             "collisionObjects will be ignored",
                             name_, o.first);
        continue;
      }
      convexes_[o.first] = {o.second.first, S_ObjectPtr(o.second.second->clone())};
      auto it = module_.collisionTransforms().find(o.first);
      if(it != module_.collisionTransforms().end())
      {
        collisionTransforms_[o.first] = it->second;
      }
      else
      {
        collisionTransforms_[o.first] = sva::PTransformd::Identity();
      }
    }
    for(const auto & b : mb().bodies())
    {
      collisionTransforms_[b.name()] = sva::PTransformd::Identity();
    }
    for(const auto & p : module_.collisionTransforms())
    {
      collisionTransforms_[p.first] = p.second;
    }
    fixCollisionTransforms();
  }

  if(loadFiles)
  {
    if(bfs::exists(module_.rsdf_dir))
    {
      loadRSDFFromDir(module_.rsdf_dir);
    }
    else if(module_.rsdf_dir.size())
    {
      mc_rtc::log::error("RSDF directory ({}) specified by RobotModule for {} does not exist.", module_.rsdf_dir,
                         module_.name);
    }
  }

  forceSensors_ = module_.forceSensors();
  if(loadFiles)
  {
    for(auto & fs : forceSensors_)
    {
      bfs::path calib_file = bfs::path(module_.calib_dir) / std::string("calib_data." + fs.name());
      fs.loadCalibrator(calib_file.string(), mbc().gravity);
    }
  }
  for(size_t i = 0; i < forceSensors_.size(); ++i)
  {
    const auto & fs = forceSensors_[i];
    forceSensorsIndex_[fs.name()] = i;
    bodyForceSensors_[fs.parentBody()] = i;
  }

  stance_ = module_.stance();

  bodySensors_ = module_.bodySensors();
  // Add a single default sensor if no sensor on the robot
  if(bodySensors_.size() == 0)
  {
    bodySensors_.emplace_back();
  }
  for(size_t i = 0; i < bodySensors_.size(); ++i)
  {
    const auto & bS = bodySensors_[i];
    bodySensorsIndex_[bS.name()] = i;
    bodyBodySensors_[bS.parentBody()] = i;
  }

  devices_ = module_.devices();
  for(size_t i = 0; i < devices_.size(); ++i)
  {
    auto & d = devices_[i];
    if(d->parent() == "")
    {
      d->parent(mb().body(0).name());
    }
    devicesIndex_[d->name()] = i;
  }

  refJointOrder_ = module_.ref_joint_order();
  refJointIndexToMBCIndex_.resize(refJointOrder_.size());
  for(size_t i = 0; i < refJointOrder_.size(); ++i)
  {
    const auto & jN = refJointOrder_[i];
    if(hasJoint(jN))
    {
      auto jIndex = mb().jointIndexByName(jN);
      refJointIndexToMBCIndex_[i] = mb().joint(jIndex).dof() != 0 ? jIndex : -1;
    }
    else
    {
      refJointIndexToMBCIndex_[i] = -1;
    }
  }

  springs_ = module_.springs();
  flexibility_ = module_.flexibility();

  zmp_ = Eigen::Vector3d::Zero();

  std::string urdf;
  auto loadUrdf = [&module_, &urdf]() -> const std::string & {
    if(urdf.size())
    {
      return urdf;
    }
    const auto & urdfPath = module_.urdf_path;
    std::ifstream ifs(urdfPath);
    if(ifs.is_open())
    {
      std::stringstream urdfSS;
      urdfSS << ifs.rdbuf();
      urdf = urdfSS.str();
      return urdf;
    }
    mc_rtc::log::error("Could not open urdf file {} for robot {}, cannot initialize grippers", urdfPath, module_.name);
    mc_rtc::log::error_and_throw("Failed to initialize grippers");
  };
  for(const auto & gripper : module_.grippers())
  {
    auto mimics = gripper.mimics();
    auto safety = gripper.safety();
    if(mimics)
    {
      grippers_[gripper.name].reset(new mc_control::Gripper(*this, gripper.joints, *mimics, gripper.reverse_limits,
                                                            safety ? *safety : module_.gripperSafety()));
    }
    else
    {
      grippers_[gripper.name].reset(new mc_control::Gripper(*this, gripper.joints, loadUrdf(), gripper.reverse_limits,
                                                            safety ? *safety : module_.gripperSafety()));
    }
  }
  for(auto & g : grippers_)
  {
    grippersRef_.push_back(std::ref(*g.second));
  }
}

const std::string & Robot::name() const
{
  return name_;
}

void Robot::name(const std::string & name)
{
  name_ = name;
}

const RobotModule & Robot::module() const
{
  return robots_->robotModule(robots_idx_);
}

BodySensor & Robot::bodySensor()
{
  return bodySensors_[0];
}

const BodySensor & Robot::bodySensor() const
{
  return bodySensors_[0];
}

bool Robot::hasBodySensor(const std::string & name) const
{
  return bodySensorsIndex_.count(name) != 0;
}

bool Robot::bodyHasBodySensor(const std::string & body) const
{
  return bodyBodySensors_.count(body) != 0;
}

BodySensor & Robot::bodySensor(const std::string & name)
{
  return const_cast<BodySensor &>(static_cast<const Robot *>(this)->bodySensor(name));
}

const BodySensor & Robot::bodySensor(const std::string & name) const
{
  return bodySensors_[bodySensorsIndex_.at(name)];
}

BodySensor & Robot::bodyBodySensor(const std::string & body)
{
  return const_cast<BodySensor &>(static_cast<const Robot *>(this)->bodyBodySensor(body));
}

const BodySensor & Robot::bodyBodySensor(const std::string & body) const
{
  return bodySensors_[bodyBodySensors_.at(body)];
}

BodySensorVector & Robot::bodySensors()
{
  return bodySensors_;
}

const BodySensorVector & Robot::bodySensors() const
{
  return bodySensors_;
}

bool Robot::hasJoint(const std::string & name) const
{
  return mb().jointIndexByName().count(name) != 0;
}

bool Robot::hasBody(const std::string & name) const
{
  return mb().bodyIndexByName().count(name) != 0;
}

unsigned int Robot::jointIndexByName(const std::string & name) const
{
  return mb().jointIndexByName().at(name);
}

int Robot::jointIndexInMBC(size_t jointIndex) const
{
  return refJointIndexToMBCIndex_.at(jointIndex);
}

unsigned int Robot::bodyIndexByName(const std::string & name) const
{
  return mb().bodyIndexByName().at(name);
}

rbd::MultiBody & Robot::mb()
{
  return robots_->mbs_[robots_idx_];
}
const rbd::MultiBody & Robot::mb() const
{
  return robots_->mbs_[robots_idx_];
}

rbd::MultiBodyConfig & Robot::mbc()
{
  return robots_->mbcs_[robots_idx_];
}
const rbd::MultiBodyConfig & Robot::mbc() const
{
  return robots_->mbcs_[robots_idx_];
}

rbd::MultiBodyGraph & Robot::mbg()
{
  return robots_->mbgs_[robots_idx_];
}
const rbd::MultiBodyGraph & Robot::mbg() const
{
  return robots_->mbgs_[robots_idx_];
}

const std::vector<std::vector<double>> & Robot::q() const
{
  return mbc().q;
}
const std::vector<std::vector<double>> & Robot::alpha() const
{
  return mbc().alpha;
}
const std::vector<std::vector<double>> & Robot::alphaD() const
{
  return mbc().alphaD;
}
const std::vector<std::vector<double>> & Robot::jointTorque() const
{
  return mbc().jointTorque;
}
const std::vector<sva::PTransformd> & Robot::bodyPosW() const
{
  return mbc().bodyPosW;
}
const std::vector<sva::MotionVecd> & Robot::bodyVelW() const
{
  return mbc().bodyVelW;
}
const std::vector<sva::MotionVecd> & Robot::bodyVelB() const
{
  return mbc().bodyVelB;
}
const std::vector<sva::MotionVecd> & Robot::bodyAccB() const
{
  return mbc().bodyAccB;
}
std::vector<std::vector<double>> & Robot::q()
{
  return mbc().q;
}
std::vector<std::vector<double>> & Robot::alpha()
{
  return mbc().alpha;
}
std::vector<std::vector<double>> & Robot::alphaD()
{
  return mbc().alphaD;
}
std::vector<std::vector<double>> & Robot::jointTorque()
{
  return mbc().jointTorque;
}
std::vector<sva::PTransformd> & Robot::bodyPosW()
{
  return mbc().bodyPosW;
}
std::vector<sva::MotionVecd> & Robot::bodyVelW()
{
  return mbc().bodyVelW;
}
std::vector<sva::MotionVecd> & Robot::bodyVelB()
{
  return mbc().bodyVelB;
}
std::vector<sva::MotionVecd> & Robot::bodyAccB()
{
  return mbc().bodyAccB;
}

const sva::PTransformd & Robot::bodyPosW(const std::string & name) const
{
  return bodyPosW()[bodyIndexByName(name)];
}

sva::PTransformd Robot::X_b1_b2(const std::string & b1, const std::string & b2) const
{
  return bodyPosW()[bodyIndexByName(b2)] * bodyPosW()[bodyIndexByName(b1)].inv();
}

const sva::MotionVecd & Robot::bodyVelW(const std::string & name) const
{
  return bodyVelW()[bodyIndexByName(name)];
}

const sva::MotionVecd & Robot::bodyVelB(const std::string & name) const
{
  return bodyVelB()[bodyIndexByName(name)];
}

const sva::MotionVecd & Robot::bodyAccB(const std::string & name) const
{
  return bodyAccB()[bodyIndexByName(name)];
}

Eigen::Vector3d Robot::com() const
{
  return rbd::computeCoM(mb(), mbc());
}
Eigen::Vector3d Robot::comVelocity() const
{
  return rbd::computeCoMVelocity(mb(), mbc());
}
Eigen::Vector3d Robot::comAcceleration() const
{
  return rbd::computeCoMAcceleration(mb(), mbc());
}

sva::ForceVecd Robot::bodyWrench(const std::string & bodyName) const
{
  if(bodyHasForceSensor(bodyName))
  { // Faster computation when there is a force sensor directly attached to the body
    const auto & fs = bodyForceSensor(bodyName);
    sva::ForceVecd w_fsactual = fs.wrenchWithoutGravity(*this);
    sva::PTransformd X_fsactual_body = fs.X_fsactual_parent();
    return X_fsactual_body.dualMul(w_fsactual);
  }
  else
  { /* If a force sensor is not directly attached to the body,
       attempt to find it up the kinematic tree */
    const auto & fs = indirectBodyForceSensor(bodyName);
    sva::ForceVecd w_fsactual = fs.wrenchWithoutGravity(*this);
    const auto & X_0_body = bodyPosW(bodyName);
    const auto & X_0_parent = bodyPosW(fs.parentBody());
    const auto X_parent_body = X_0_body * X_0_parent.inv();
    sva::PTransformd X_fsactual_body = X_parent_body * fs.X_fsactual_parent();
    return X_fsactual_body.dualMul(w_fsactual);
  }
}

sva::ForceVecd Robot::surfaceWrench(const std::string & surfaceName) const
{
  const auto & surface = this->surface(surfaceName);
  const auto & bodyName = surface.bodyName();
  return surface.X_b_s().dualMul(bodyWrench(bodyName));
}

Eigen::Vector2d Robot::cop(const std::string & surfaceName, double min_pressure) const
{
  const sva::ForceVecd w_surf = surfaceWrench(surfaceName);
  const double pressure = w_surf.force()(2);
  if(pressure < min_pressure)
  {
    return Eigen::Vector2d::Zero();
  }
  const Eigen::Vector3d & tau_surf = w_surf.couple();
  return Eigen::Vector2d(-tau_surf(1) / pressure, +tau_surf(0) / pressure);
}

Eigen::Vector3d Robot::copW(const std::string & surfaceName, double min_pressure) const
{
  Eigen::Vector3d cop_s;
  cop_s << cop(surfaceName, min_pressure), 0.;
  const sva::PTransformd X_0_s = surfacePose(surfaceName);
  return X_0_s.translation() + X_0_s.rotation().transpose() * cop_s;
}

sva::ForceVecd Robot::netWrench(const std::vector<std::string> & sensorNames) const
{
  // Compute net total wrench from all sensors in contact
  sva::ForceVecd netTotalWrench{sva::ForceVecd::Zero()};
  for(const auto & sensorName : sensorNames)
  {
    const auto & sensor = forceSensor(sensorName);
    netTotalWrench += sensor.worldWrenchWithoutGravity(*this);
  }
  return netTotalWrench;
}

Eigen::Vector3d Robot::zmp(const sva::ForceVecd & netTotalWrench,
                           const Eigen::Vector3d & plane_p,
                           const Eigen::Vector3d & plane_n,
                           double minimalNetNormalForce) const
{
  return mc_rbdyn::zmp(netTotalWrench, plane_p, plane_n, minimalNetNormalForce);
}

Eigen::Vector3d Robot::zmp(const sva::ForceVecd & netWrench,
                           const sva::PTransformd & zmpFrame,
                           double minimalNetNormalForce) const
{
  return mc_rbdyn::zmp(netWrench, zmpFrame, minimalNetNormalForce);
}

bool Robot::maybeZMP(const sva::ForceVecd & netWrench,
                     const sva::PTransformd & zmpFrame,
                     Eigen::Vector3d & result,
                     double minimalNetNormalForce) const
{
  return mc_rbdyn::maybeZMP(netWrench, zmpFrame, result, minimalNetNormalForce);
}

Eigen::Vector3d Robot::zmp(const std::vector<std::string> & sensorNames,
                           const Eigen::Vector3d & plane_p,
                           const Eigen::Vector3d & plane_n,
                           double minimalNetNormalForce) const
{
  return zmp(netWrench(sensorNames), plane_p, plane_n, minimalNetNormalForce);
}

Eigen::Vector3d Robot::zmp(const std::vector<std::string> & sensorNames,
                           const sva::PTransformd & zmpFrame,
                           double minimalNetNormalForce) const
{
  Eigen::Vector3d n = zmpFrame.rotation().row(2);
  Eigen::Vector3d p = zmpFrame.translation();
  return zmp(sensorNames, p, n, minimalNetNormalForce);
}

const std::vector<std::vector<double>> & Robot::ql() const
{
  return ql_;
}
const std::vector<std::vector<double>> & Robot::qu() const
{
  return qu_;
}
const std::vector<std::vector<double>> & Robot::vl() const
{
  return vl_;
}
const std::vector<std::vector<double>> & Robot::vu() const
{
  return vu_;
}
const std::vector<std::vector<double>> & Robot::al() const
{
  return al_;
}
const std::vector<std::vector<double>> & Robot::au() const
{
  return au_;
}
const std::vector<std::vector<double>> & Robot::jl() const
{
  return jl_;
}
const std::vector<std::vector<double>> & Robot::ju() const
{
  return ju_;
}
const std::vector<std::vector<double>> & Robot::tl() const
{
  return tl_;
}
const std::vector<std::vector<double>> & Robot::tu() const
{
  return tu_;
}
const std::vector<std::vector<double>> & Robot::tdl() const
{
  return tdl_;
}
const std::vector<std::vector<double>> & Robot::tdu() const
{
  return tdu_;
}
std::vector<std::vector<double>> & Robot::ql()
{
  return ql_;
}
std::vector<std::vector<double>> & Robot::qu()
{
  return qu_;
}
std::vector<std::vector<double>> & Robot::vl()
{
  return vl_;
}
std::vector<std::vector<double>> & Robot::vu()
{
  return vu_;
}
std::vector<std::vector<double>> & Robot::al()
{
  return al_;
}
std::vector<std::vector<double>> & Robot::au()
{
  return au_;
}
std::vector<std::vector<double>> & Robot::jl()
{
  return jl_;
}
std::vector<std::vector<double>> & Robot::ju()
{
  return ju_;
}
std::vector<std::vector<double>> & Robot::tl()
{
  return tl_;
}
std::vector<std::vector<double>> & Robot::tu()
{
  return tu_;
}
std::vector<std::vector<double>> & Robot::tdl()
{
  return tdl_;
}
std::vector<std::vector<double>> & Robot::tdu()
{
  return tdu_;
}

const std::vector<Flexibility> & Robot::flexibility() const
{
  return flexibility_;
}

std::vector<Flexibility> & Robot::flexibility()
{
  return flexibility_;
}

const std::vector<double> & Robot::encoderValues() const
{
  return encoderValues_;
}

void Robot::encoderValues(const std::vector<double> & encoderValues)
{
  encoderValues_ = encoderValues;
}

const std::vector<double> & Robot::encoderVelocities() const
{
  return encoderVelocities_;
}

void Robot::encoderVelocities(const std::vector<double> & encoderVelocities)
{
  encoderVelocities_ = encoderVelocities;
}

const std::vector<double> & Robot::flexibilityValues() const
{
  return flexibilityValues_;
}

void Robot::flexibilityValues(const std::vector<double> & flexibilityValues)
{
  flexibilityValues_ = flexibilityValues;
}

const std::vector<double> & Robot::jointTorques() const
{
  return jointTorques_;
}

void Robot::jointTorques(const std::vector<double> & jointTorques)
{
  jointTorques_ = jointTorques;
}

const std::vector<std::string> & Robot::refJointOrder() const
{
  return refJointOrder_;
}

bool Robot::hasForceSensor(const std::string & name) const
{
  return forceSensorsIndex_.count(name) != 0;
}

bool Robot::bodyHasForceSensor(const std::string & body) const
{
  return bodyForceSensors_.count(body) != 0;
}

bool Robot::bodyHasIndirectForceSensor(const std::string & body) const
{
  return bodyHasForceSensor(body) || findIndirectForceSensorBodyName(body).size();
}

bool Robot::surfaceHasForceSensor(const std::string & surfaceName) const
{
  return bodyHasForceSensor(surface(surfaceName).bodyName());
}

bool Robot::surfaceHasIndirectForceSensor(const std::string & surfaceName) const
{
  return bodyHasIndirectForceSensor(surface(surfaceName).bodyName());
}

ForceSensor & Robot::forceSensor(const std::string & name)
{
  return const_cast<ForceSensor &>(static_cast<const Robot *>(this)->forceSensor(name));
}

const ForceSensor & Robot::forceSensor(const std::string & name) const
{
  return forceSensors_[forceSensorsIndex_.at(name)];
}

ForceSensor & Robot::bodyForceSensor(const std::string & body)
{
  return const_cast<ForceSensor &>(static_cast<const Robot *>(this)->bodyForceSensor(body));
}

const ForceSensor & Robot::bodyForceSensor(const std::string & body) const
{
  return forceSensors_.at(bodyForceSensors_.at(body));
}

ForceSensor & Robot::surfaceForceSensor(const std::string & surfaceName)
{
  return bodyForceSensor(surface(surfaceName).bodyName());
}

const ForceSensor & Robot::surfaceForceSensor(const std::string & surfaceName) const
{
  return bodyForceSensor(surface(surfaceName).bodyName());
}

std::string Robot::findIndirectForceSensorBodyName(const std::string & body) const
{
  int nextIndex = mb().bodyIndexByName().at(body);
  while(nextIndex >= 0)
  {
    const auto & b = mb().body(nextIndex);
    if(bodyHasForceSensor(b.name()))
    {
      return b.name();
    }
    nextIndex = mb().parent(nextIndex);
  }
  return std::string{};
}

const ForceSensor & Robot::indirectBodyForceSensor(const std::string & body) const
{
  const auto bodyName = findIndirectForceSensorBodyName(body);
  if(bodyName.empty())
  {
    mc_rtc::log::error_and_throw("No force sensor (directly or indirectly) attached to body {}", body);
  }
  return bodyForceSensor(bodyName);
}

ForceSensor & Robot::indirectBodyForceSensor(const std::string & body)
{
  return const_cast<ForceSensor &>(static_cast<const Robot *>(this)->indirectBodyForceSensor(body));
}

const ForceSensor & Robot::indirectSurfaceForceSensor(const std::string & surfaceName) const
{
  return indirectBodyForceSensor(surface(surfaceName).bodyName());
}

ForceSensor & Robot::indirectSurfaceForceSensor(const std::string & surface)
{
  return const_cast<ForceSensor &>(static_cast<const Robot *>(this)->indirectSurfaceForceSensor(surface));
}

bool Robot::hasSurface(const std::string & surface) const
{
  return surfaces_.count(surface) != 0;
}

std::vector<ForceSensor> & Robot::forceSensors()
{
  return forceSensors_;
}

const std::vector<ForceSensor> & Robot::forceSensors() const
{
  return forceSensors_;
}

mc_rbdyn::Surface & Robot::surface(const std::string & sName)
{
  return const_cast<mc_rbdyn::Surface &>(static_cast<const Robot *>(this)->surface(sName));
}

sva::PTransformd Robot::surfacePose(const std::string & sName) const
{
  return surface(sName).X_0_s(*this);
}

const mc_rbdyn::Surface & Robot::surface(const std::string & sName) const
{
  if(!hasSurface(sName))
  {
    mc_rtc::log::error_and_throw("No surface named {} found in robot {}", sName, this->name());
  }
  return *(surfaces_.at(sName));
}

const std::map<std::string, SurfacePtr> & Robot::surfaces() const
{
  return surfaces_;
}

std::vector<std::string> Robot::availableSurfaces() const
{
  std::vector<std::string> ret;
  ret.reserve(surfaces_.size());
  for(const auto & s : surfaces_)
  {
    ret.push_back(s.first);
  }
  return ret;
}

bool Robot::hasConvex(const std::string & name) const
{
  return convexes_.count(name);
}

Robot::convex_pair_t & Robot::convex(const std::string & cName)
{
  return const_cast<Robot::convex_pair_t &>(static_cast<const Robot *>(this)->convex(cName));
}
const Robot::convex_pair_t & Robot::convex(const std::string & cName) const
{
  if(convexes_.count(cName) == 0)
  {
    mc_rtc::log::error_and_throw("No convex named {} found in robot {}", cName, this->name_);
  }
  return convexes_.at(cName);
}

const std::map<std::string, Robot::convex_pair_t> & Robot::convexes() const
{
  return convexes_;
}

void Robot::addConvex(const std::string & cName,
                      const std::string & body,
                      S_ObjectPtr convex,
                      const sva::PTransformd & X_b_c)
{
  if(convexes_.count(cName))
  {
    mc_rtc::log::error("Attempted to add a convex named {} that already exists in {}", cName, name());
    return;
  }
  convexes_[cName] = {body, convex};
  collisionTransforms_[cName] = X_b_c;
  sch::mc_rbdyn::transform(*convex, X_b_c * bodyPosW(body));
}

void Robot::removeConvex(const std::string & cName)
{
  if(convexes_.count(cName))
  {
    convexes_.erase(cName);
    collisionTransforms_.erase(cName);
  }
}

const sva::PTransformd & Robot::bodyTransform(const std::string & bName) const
{
  if(!hasBody(bName))
  {
    mc_rtc::log::error_and_throw("No body transform with name {} found in this robot", bName);
  }
  return bodyTransforms_[bodyIndexByName(bName)];
}

const sva::PTransformd & Robot::bodyTransform(int bodyIndex) const
{
  return bodyTransforms_[bodyIndex];
}

const std::vector<sva::PTransformd> & Robot::bodyTransforms() const
{
  return bodyTransforms_;
}

const sva::PTransformd & Robot::collisionTransform(const std::string & cName) const
{
  if(collisionTransforms_.count(cName) == 0)
  {
    mc_rtc::log::error_and_throw("No collision transform with name {} found in this robot", cName);
  }
  return collisionTransforms_.at(cName);
}

void Robot::fixSurfaces()
{
  for(auto & s : surfaces_)
  {
    const sva::PTransformd & trans = bodyTransform(s.second->bodyName());
    s.second->X_b_s(s.second->X_b_s() * trans);
  }
}

void Robot::fixCollisionTransforms()
{
  for(auto & ct : collisionTransforms_)
  {
    if(convexes_.count(ct.first))
    {
      const auto & trans = bodyTransform(convexes_.at(ct.first).first);
      ct.second = ct.second * trans;
    }
    else
    {
      const auto & trans = bodyTransform(ct.first);
      ct.second = ct.second * trans;
    }
  }
}

void Robot::loadRSDFFromDir(const std::string & surfaceDir)
{
  std::vector<SurfacePtr> surfacesIn = readRSDFFromDir(surfaceDir);
  for(const auto & sp : surfacesIn)
  {
    /* Check coherence of surface with mb */
    if(hasBody(sp->bodyName()))
    {
      surfaces_[sp->name()] = sp;
    }
    else
    {
      mc_rtc::log::warning("Loaded surface {} attached to body {} from RSDF but the robot {} has no such body, discard "
                           "this surface to avoid future problems...",
                           sp->name(), sp->bodyName(), name());
    }
  }
  fixSurfaces();
}

std::map<std::string, std::vector<double>> Robot::stance() const
{
  return stance_;
}

unsigned int mc_rbdyn::Robot::robotIndex() const
{
  return robots_idx_;
}

void Robot::forwardKinematics()
{
  rbd::forwardKinematics(mb(), mbc());
}
void Robot::forwardKinematics(rbd::MultiBodyConfig & mbc) const
{
  rbd::forwardKinematics(mb(), mbc);
}

void Robot::forwardVelocity()
{
  rbd::forwardVelocity(mb(), mbc());
}
void Robot::forwardVelocity(rbd::MultiBodyConfig & mbc) const
{
  rbd::forwardVelocity(mb(), mbc);
}

void Robot::forwardAcceleration(const sva::MotionVecd & A_0)
{
  rbd::forwardAcceleration(mb(), mbc(), A_0);
}
void Robot::forwardAcceleration(rbd::MultiBodyConfig & mbc, const sva::MotionVecd & A_0) const
{
  rbd::forwardAcceleration(mb(), mbc, A_0);
}

void mc_rbdyn::Robot::eulerIntegration(double step)
{
  rbd::eulerIntegration(mb(), mbc(), step);
}

void mc_rbdyn::Robot::eulerIntegration(rbd::MultiBodyConfig & mbc, double step) const
{
  rbd::eulerIntegration(mb(), mbc, step);
}

const sva::PTransformd & Robot::posW() const
{
  return bodyPosW().at(0);
}

void Robot::posW(const sva::PTransformd & pt)
{
  if(mb().joint(0).type() == rbd::Joint::Type::Free)
  {
    Eigen::Quaterniond rotation{pt.rotation().transpose()};
    rotation.normalize();
    q()[0] = {rotation.w(),         rotation.x(),         rotation.y(),        rotation.z(),
              pt.translation().x(), pt.translation().y(), pt.translation().z()};
    forwardKinematics();
  }
  else if(mb().joint(0).type() == rbd::Joint::Type::Fixed)
  {
    sva::PTransformd pt_ = pt;
    pt_.rotation() = Eigen::Quaterniond(pt.rotation()).normalized().toRotationMatrix();
    mb().transform(0, pt_);
    forwardKinematics();
    fixSCH(*this, this->convexes_, this->collisionTransforms_);
  }
  else
  {
    mc_rtc::log::error_and_throw<std::logic_error>(
        "The root pose can only be changed for robots with a free flyer or a fixed joint as joint(0)");
  }
}

void Robot::velW(const sva::MotionVecd & vel)
{
  if(mb().joint(0).type() == rbd::Joint::Type::Free)
  {
    auto vB = sva::PTransformd(mbc().bodyPosW[0].rotation()) * vel;
    alpha()[0][0] = vB.angular().x();
    alpha()[0][1] = vB.angular().y();
    alpha()[0][2] = vB.angular().z();
    alpha()[0][3] = vB.linear().x();
    alpha()[0][4] = vB.linear().y();
    alpha()[0][5] = vB.linear().z();
    forwardVelocity();
  }
  else
  {
    mc_rtc::log::warning("You cannot set the base velocity on a fixed-base robot");
  }
}

const sva::MotionVecd & Robot::velW() const
{
  return bodyVelW().at(0);
}

void Robot::accW(const sva::MotionVecd & acc)
{
  if(mb().joint(0).type() == rbd::Joint::Type::Free)
  {
    auto aB = sva::PTransformd(mbc().bodyPosW[0].rotation()) * acc;
    alphaD()[0][0] = aB.angular().x();
    alphaD()[0][1] = aB.angular().y();
    alphaD()[0][2] = aB.angular().z();
    alphaD()[0][3] = aB.linear().x();
    alphaD()[0][4] = aB.linear().y();
    alphaD()[0][5] = aB.linear().z();
    forwardAcceleration();
  }
  else
  {
    mc_rtc::log::warning("You cannot set the base acceleration on a fixed-base robot");
  }
}

const sva::MotionVecd Robot::accW() const
{
  Eigen::Matrix3d rot = posW().rotation().transpose();
  return sva::PTransformd{rot} * mbc().bodyAccB[0];
}

void Robot::copyLoadedData(Robot & robot) const
{
  for(const auto & s : surfaces_)
  {
    robot.surfaces_[s.first] = s.second->copy();
  }
  robot.fixSurfaces();
  for(const auto & cH : convexes_)
  {
    robot.convexes_[cH.first] = {cH.second.first, S_ObjectPtr(cH.second.second->clone())};
  }
  robot.collisionTransforms_ = collisionTransforms_;
  robot.fixCollisionTransforms();
  fixSCH(robot, robot.convexes_, robot.collisionTransforms_);
  for(size_t i = 0; i < forceSensors_.size(); ++i)
  {
    robot.forceSensors_[i].copyCalibrator(forceSensors_[i]);
  }
}

mc_rbdyn::Surface & Robot::copySurface(const std::string & sName, const std::string & name)
{
  if(hasSurface(name))
  {
    mc_rtc::log::error_and_throw("{} already exists within robot {}. Cannot overwrite an existing surface", name,
                                 this->name_);
  }
  const Surface & surf = surface(sName);
  SurfacePtr nSurf = surf.copy();
  nSurf->name(name);
  surfaces_[name] = nSurf;
  return *nSurf;
}

void mc_rbdyn::Robot::addSurface(SurfacePtr surface, bool doNotReplace)
{
  if(!hasBody(surface->bodyName()))
  {
    mc_rtc::log::warning("Surface {} attached to body {} but robot {} has no such body.", surface->name(),
                         surface->bodyName(), name());
    return;
  }
  if(hasSurface(surface->name()) && doNotReplace)
  {
    mc_rtc::log::warning("Surface {} already exists for the robot {}.", surface->name(), name());
    return;
  }
  surfaces_[surface->name()] = std::move(surface);
}

MC_RTC_diagnostic_pop

double mc_rbdyn::Robot::mass() const
{
  double mass = 0.;
  for(const auto & b : mb().bodies())
  {
    mass += b.inertia().mass();
  }
  return mass;
}

void mc_rbdyn::Robot::zmpTarget(const Eigen::Vector3d & zmp)
{
  zmp_ = zmp;
}

const Eigen::Vector3d & mc_rbdyn::Robot::zmpTarget() const
{
  return zmp_;
}

mc_control::Gripper & Robot::gripper(const std::string & gripper)
{
  if(!grippers_.count(gripper))
  {
    mc_rtc::log::error_and_throw("No gripper named {} in robot {}", gripper, name());
  }
  return *grippers_.at(gripper);
}

bool Robot::hasGripper(const std::string & gripper) const
{
  return grippers_.count(gripper);
}

unsigned int robotIndexFromConfig(const mc_rtc::Configuration & config,
                                  const mc_rbdyn::Robots & robots,
                                  const std::string & prefix,
                                  bool required,
                                  const std::string & robotIndexKey,
                                  const std::string & robotNameKey,
                                  const std::string & defaultRobotName)
{
  const auto & robot = robotFromConfig(config, robots, prefix, required, robotIndexKey, robotNameKey, defaultRobotName);
  return robot.robotIndex();
}

std::string robotNameFromConfig(const mc_rtc::Configuration & config,
                                const mc_rbdyn::Robots & robots,
                                const std::string & prefix,
                                bool required,
                                const std::string & robotIndexKey,
                                const std::string & robotNameKey,
                                const std::string & defaultRobotName)
{
  const auto & robot = robotFromConfig(config, robots, prefix, required, robotIndexKey, robotNameKey, defaultRobotName);
  return robot.name();
}

const mc_rbdyn::Robot & robotFromConfig(const mc_rtc::Configuration & config,
                                        const mc_rbdyn::Robots & robots,
                                        const std::string & prefix,
                                        bool required,
                                        const std::string & robotIndexKey,
                                        const std::string & robotNameKey,
                                        const std::string & defaultRobotName)
{
  auto p = std::string{""};
  if(prefix.size())
  {
    p = "[" + prefix + "] ";
  }
  if(config.has(robotNameKey))
  {
    const std::string & robotName = config(robotNameKey);
    if(robots.hasRobot(robotName))
    {
      return robots.robot(robotName);
    }
    else
    {
      mc_rtc::log::error_and_throw("{} No robot named {} in this controller", p, robotName);
    }
  }
  else if(config.has(robotIndexKey))
  {
    mc_rtc::log::warning("[MC_RTC_DEPRECATED]{} \"robotIndex\" will be deprecated in future versions, use \"robot: "
                         "<robot name>\" instead",
                         p);
    const unsigned int robotIndex = config(robotIndexKey);
    if(robotIndex < robots.size())
    {
      return robots.robot(robotIndex);
    }
    else
    {
      mc_rtc::log::error_and_throw("{}No robot with index {} in this controller ({} robots loaded)", p, robotIndex,
                                   robots.size());
    }
  }
  else
  {
    if(!required)
    {
      return defaultRobotName.size() ? robots.robot(defaultRobotName) : robots.robot();
    }
    else
    {
      mc_rtc::log::error_and_throw("{} \"robotName\" is required.", p);
    }
  }
}

void Robot::addDevice(DevicePtr device)
{
  if(devicesIndex_.count(device->name()))
  {
    mc_rtc::log::error_and_throw("You cannot have multiple generic sensor with the same name in a robot");
  }
  devices_.push_back(std::move(device));
  auto & d = devices_.back();
  if(d->parent() == "")
  {
    d->parent(mb().body(0).name());
  }
  devicesIndex_[device->name()] = devices_.size() - 1;
}

} // namespace mc_rbdyn

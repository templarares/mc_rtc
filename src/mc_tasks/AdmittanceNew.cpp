/*
 * Copyright 2015-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_tasks/AdmittanceNew.h>

#include <mc_tasks/MetaTaskLoader.h>

#include <mc_rbdyn/configuration_io.h>
#include <mc_rbdyn/rpy_utils.h>

#include <mc_rtc/gui/ArrayLabel.h>
#include <mc_rtc/gui/Transform.h>

namespace mc_tasks
{

namespace force
{

using mc_filter::utils::clampInPlaceAndWarn;

AdmittanceNew::AdmittanceNew(const std::string & surfaceName,
                               const mc_rbdyn::Robots & robots,
                               unsigned int robotIndex,
                               double stiffness,
                               double weight)
: SurfaceTransformTask(surfaceName, robots, robotIndex, stiffness, weight), robots_(robots), rIndex_(robotIndex),
  surface_(robots.robot(robotIndex).surface(surfaceName))
{
  const auto & robot = robots.robot(robotIndex);
  if(!robot.surfaceHasIndirectForceSensor(surfaceName))
  {
    // mc_rtc::log::error_and_throw<std::runtime_error>(
    //     "[mc_tasks::AdmittanceNew] Surface {} does not have a force sensor attached", surfaceName);
    mc_rtc::log::info("No force sensor attached. Use other force sensors");
    mc_rtc::log::info("Other force sensors are: {},{},{}", robot.forceSensors()[0].name(),robot.forceSensors()[1].name(),robot.forceSensors()[2].name());
    this->deducedForce(true);
  }
  name_ = "admittance_" + robots_.robot(rIndex_).name() + "_" + surfaceName;
  reset();
}

void AdmittanceNew::update(mc_solver::QPSolver &)
{
  // Compute wrench error
  wrenchError_ = measuredWrench() - targetWrench_;

  // Compute linear and angular velocity based on wrench error and admittance
  Eigen::Vector3d linearVel = admittance_.force().cwiseProduct(wrenchError_.force());
  Eigen::Vector3d angularVel = admittance_.couple().cwiseProduct(wrenchError_.couple());

  // Clamp both values in order to have a 'security'
  clampInPlaceAndWarn(linearVel, (-maxLinearVel_).eval(), maxLinearVel_, name_ + " linear velocity");
  clampInPlaceAndWarn(angularVel, (-maxAngularVel_).eval(), maxAngularVel_, name_ + " angular velocity");

  // Filter
  refVelB_ = velFilterGain_ * refVelB_ + (1 - velFilterGain_) * sva::MotionVecd(angularVel, linearVel);

  // Compute position and rotation delta
  sva::PTransformd delta(mc_rbdyn::rpyToMat(timestep_ * refVelB_.angular()), timestep_ * refVelB_.linear());

  // Apply feed forward term
  refVelB_ += feedforwardVelB_;

  // Acceleration
  SurfaceTransformTask::refAccel((refVelB_ - SurfaceTransformTask::refVelB()) / timestep_);

  // Velocity
  SurfaceTransformTask::refVelB(refVelB_);

  // Position
  target(delta * target());
}

std::function<bool(const mc_tasks::MetaTask &, std::string &)> AdmittanceNew::buildCompletionCriteria(
    double dt,
    const mc_rtc::Configuration & config) const
{
  if(config.has("wrench"))
  {
    if(robots.robot(rIndex).surfaceHasIndirectForceSensor(surfaceName))
    {
      // mc_rtc::log::error_and_throw<std::invalid_argument>("[{}] Attempted to use \"wrench\" as completion criteria but "
      //                                                     "surface \"{}\" is not attached to a force sensor",
      //                                                     name(), surfaceName);
      
      return SurfaceTransformTask::buildCompletionCriteria(dt, config);
    }
    sva::ForceVecd target_w = config("wrench");
    Eigen::Vector6d target = target_w.vector();
    Eigen::Vector6d dof = Eigen::Vector6d::Ones();
    for(int i = 0; i < 6; ++i)
    {
      if(std::isnan(target(i)))
      {
        dof(i) = 0.;
        target(i) = 0.;
      }
      else if(target(i) < 0)
      {
        dof(i) = -1.;
      }
    }
    return [this, dof, target](const mc_tasks::MetaTask & t, std::string & out) {
      const auto & self = static_cast<const mc_tasks::SurfaceTransformTask &>(t);
      //Eigen::Vector6d w = self.robots.robot(self.rIndex).surfaceWrench(self.surface()).vector();
      std::vector<std::string> contactSensors;      
      contactSensors.push_back("LeftFootForceSensor");
      contactSensors.push_back("RightFootForceSensor");
      contactSensors.push_back("LeftHandForceSensor");
      mc_rtc::log::info("checking wrench...");
      auto res=robots_.robot(rIndex_).netWrench(contactSensors);
      res*=-1.;
      res = surface_.X_0_s(robots_.robot(rIndex_)).dualMul(res);
      auto w = res.vector();
      mc_rtc::log::info("Completion check:deduced hip wrench in surface frame, assuming zero net torque, is [{},{},{},{},{},{}]", w[0],w[1],w[2],w[3],w[4],w[5]);
      for(int i = 0; i < 6; ++i)
      {
        if(dof(i) * fabs(w(i)) < target(i))
        {
          return false;
        }
      }
      out += "wrench";
      mc_rtc::log::info("wrench check passed");
      return true;
    };
  }
  return MetaTask::buildCompletionCriteria(dt, config);
}


void AdmittanceNew::reset()
{
  SurfaceTransformTask::reset();
  admittance_ = sva::ForceVecd(Eigen::Vector6d::Zero());
  feedforwardVelB_ = sva::MotionVecd(Eigen::Vector6d::Zero());
  targetWrench_ = sva::ForceVecd(Eigen::Vector6d::Zero());
  wrenchError_ = sva::ForceVecd(Eigen::Vector6d::Zero());
}

/*! \brief Load parameters from a Configuration object */
void AdmittanceNew::load(mc_solver::QPSolver & solver, const mc_rtc::Configuration & config)
{
  if(config.has("admittance"))
  {
    admittance(config("admittance"));
  }
  else if(config.has("targetPose"))
  {
    mc_rtc::log::warning("[{}] property \"targetPose\" is deprecated, use \"target\" instead", name());
    targetPose(config("targetPose"));
  }
  if(config.has("wrench"))
  {
    targetWrench(config("wrench"));
  }
  if(config.has("refVelB"))
  {
    refVelB(config("refVelB"));
  }
  if(config.has("maxVel"))
  {
    sva::MotionVecd maxVel = config("maxVel");
    maxLinearVel(maxVel.linear());
    maxAngularVel(maxVel.angular());
  }
  SurfaceTransformTask::load(solver, config);
}

void AdmittanceNew::addToLogger(mc_rtc::Logger & logger)
{
  SurfaceTransformTask::addToLogger(logger);
  MC_RTC_LOG_HELPER(name_ + "_admittance", admittance_);
  MC_RTC_LOG_HELPER(name_ + "_measured_wrench", measuredWrench);
  MC_RTC_LOG_HELPER(name_ + "_target_body_vel", feedforwardVelB_);
  MC_RTC_LOG_HELPER(name_ + "_target_wrench", targetWrench_);
  MC_RTC_LOG_HELPER(name_ + "_vel_filter_gain", velFilterGain_);
}

void AdmittanceNew::addToGUI(mc_rtc::gui::StateBuilder & gui)
{
  gui.addElement(
      {"Tasks", name_},
      mc_rtc::gui::Transform(
          "pos_target", [this]() { return this->targetPose(); },
          [this](const sva::PTransformd & pos) { this->targetPose(pos); }),
      mc_rtc::gui::Transform(
          "pos", [this]() { return robots.robot(rIndex).surface(surfaceName).X_0_s(robots.robot(rIndex)); }),
      mc_rtc::gui::ArrayInput(
          "admittanceNew", {"cx", "cy", "cz", "fx", "fy", "fz"}, [this]() { return this->admittance().vector(); },
          [this](const Eigen::Vector6d & a) { this->admittance(a); }),
      mc_rtc::gui::ArrayInput(
          "wrench", {"cx", "cy", "cz", "fx", "fy", "fz"}, [this]() { return this->targetWrench().vector(); },
          [this](const Eigen::Vector6d & a) { this->targetWrench(a); }),
      mc_rtc::gui::ArrayLabel("measured_wrench", {"cx", "cy", "cz", "fx", "fy", "fz"},
                              [this]() { return this->measuredWrench().vector(); }),
      mc_rtc::gui::NumberInput(
          "Velocity filter gain", [this]() { return velFilterGain_; }, [this](double g) { velFilterGain(g); }));
  // Don't add SurfaceTransformTask as target configuration is different
  TrajectoryTaskGeneric<tasks::qp::SurfaceTransformTask>::addToGUI(gui);
}

void AdmittanceNew::addToSolver(mc_solver::QPSolver & solver)
{
  timestep_ = solver.dt();
  SurfaceTransformTask::addToSolver(solver);
}

} // namespace force

} // namespace mc_tasks

namespace
{

static auto registered = mc_tasks::MetaTaskLoader::register_load_function(
    "admittanceNew",
    [](mc_solver::QPSolver & solver, const mc_rtc::Configuration & config) {
      auto t = std::make_shared<mc_tasks::force::AdmittanceNew>(
          config("surface"), solver.robots(), robotIndexFromConfig(config, solver.robots(), "admittanceNew"));

      t->reset();
      t->load(solver, config);
      return t;
    });
}

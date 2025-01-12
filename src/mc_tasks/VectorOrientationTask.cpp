/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_tasks/VectorOrientationTask.h>

#include <mc_tasks/MetaTaskLoader.h>

#include <mc_rtc/gui/ArrayInput.h>
#include <mc_rtc/gui/Arrow.h>
#include <mc_rtc/gui/Point3D.h>

namespace mc_tasks
{

inline static tasks::qp::VectorOrientationTask * tasks_error(mc_rtc::void_ptr & ptr)
{
  return static_cast<tasks::qp::VectorOrientationTask *>(ptr.get());
}

inline static const tasks::qp::VectorOrientationTask * tasks_error(const mc_rtc::void_ptr & ptr)
{
  return static_cast<const tasks::qp::VectorOrientationTask *>(ptr.get());
}

VectorOrientationTask::VectorOrientationTask(const std::string & bodyName,
                                             const Eigen::Vector3d & bodyVector,
                                             const Eigen::Vector3d & targetVector,
                                             const mc_rbdyn::Robots & robots,
                                             unsigned int robotIndex,
                                             double stiffness,
                                             double weight)
: VectorOrientationTask(robots.robot(robotIndex).frame(bodyName), bodyVector, stiffness, weight)
{
  this->targetVector(targetVector);
}

VectorOrientationTask::VectorOrientationTask(const mc_rbdyn::RobotFrame & frame,
                                             const Eigen::Vector3d & frameVector,
                                             double stiffness,
                                             double weight)
: TrajectoryTaskGeneric(frame, stiffness, weight), frame_(frame)
{
  const auto & X_b_f = frame.X_b_f();
  Eigen::Vector3d bodyVector = (sva::PTransformd{frameVector} * X_b_f).translation().normalized();
  switch(backend_)
  {
    case Backend::Tasks:
      finalize<tasks::qp::VectorOrientationTask>(robots.mbs(), static_cast<int>(rIndex), frame.body(), bodyVector,
                                                 bodyVector);
      break;
    default:
      mc_rtc::log::error_and_throw("[VectorOrientationTask] Not implemented for solver backend: {}", backend_);
  }
  type_ = "vectorOrientation";
  name_ = "vector_orientation_" + frame.robot().name() + "_" + frame.name();
  reset();
}

VectorOrientationTask::VectorOrientationTask(const std::string & bodyName,
                                             const Eigen::Vector3d & bodyVector,
                                             const mc_rbdyn::Robots & robots,
                                             unsigned int robotIndex,
                                             double stiffness,
                                             double weight)
: VectorOrientationTask(robots.robot(robotIndex).frame(bodyName), bodyVector, stiffness, weight)
{
}

void VectorOrientationTask::reset()
{
  TrajectoryTaskGeneric::reset();
  // Should be errorT->actual(), but it is not computed until the first call to
  // errorT::update()
  switch(backend_)
  {
    case Backend::Tasks:
    {
      Eigen::Matrix3d E_0_b = frame_->robot().frame(body()).position().rotation().transpose();
      Eigen::Vector3d actualVector = E_0_b * tasks_error(errorT)->bodyVector();
      this->targetVector(actualVector.normalized());
    }
    break;
    default:
      break;
  }
}

void VectorOrientationTask::targetVector(const Eigen::Vector3d & ori)
{
  switch(backend_)
  {
    case Backend::Tasks:
    {
      tasks_error(errorT)->target((sva::PTransformd{ori} * frame_->X_b_f()).translation().normalized());
    }
    break;
    default:
      break;
  }
}

Eigen::Vector3d VectorOrientationTask::targetVector() const
{
  switch(backend_)
  {
    case Backend::Tasks:
    {
      return (frame_->X_b_f() * tasks_error(errorT)->target()).translation();
    }
    default:
      mc_rtc::log::error_and_throw("Not implemented");
  }
}

Eigen::Vector3d VectorOrientationTask::actual() const
{
  switch(backend_)
  {
    case Backend::Tasks:
    {
      return (frame_->X_b_f() * tasks_error(errorT)->actual()).translation();
    }
    default:
      mc_rtc::log::error_and_throw("Not implemented");
  }
}

void VectorOrientationTask::addToLogger(mc_rtc::Logger & logger)
{
  TrajectoryBase::addToLogger(logger);
  MC_RTC_LOG_GETTER(name_ + "_target", targetVector);
  MC_RTC_LOG_HELPER(name_ + "_current", actual);
  MC_RTC_LOG_HELPER(name_ + "_error", eval);
}

void VectorOrientationTask::addToGUI(mc_rtc::gui::StateBuilder & gui)
{
  TrajectoryBase::addToGUI(gui);
  gui.addElement(
      {"Tasks", name_},
      mc_rtc::gui::ArrayInput(
          "Target Direction", {"x", "y", "z"}, [this]() { return targetVector(); },
          [this](const Eigen::Vector3d & target) { targetVector(target); }),
      mc_rtc::gui::Arrow(
          "Actual", mc_rtc::gui::ArrowConfig(mc_rtc::gui::Color(0., 0., 1.)),
          [this]() -> Eigen::Vector3d { return frame_->position().translation(); },
          [this]() -> Eigen::Vector3d { return frame_->position().translation() + .25 * actual(); }),
      mc_rtc::gui::Arrow(
          "Target", mc_rtc::gui::ArrowConfig(mc_rtc::gui::Color(1., 0., 0.)),
          [this]() -> Eigen::Vector3d { return frame_->position().translation(); },
          [this]() -> Eigen::Vector3d { return frame_->position().translation() + .25 * targetVector(); }),
      mc_rtc::gui::Point3D(
          "Arrow end point",
          [this]() -> Eigen::Vector3d { return frame_->position().translation() + .25 * targetVector(); },
          [this](const Eigen::Vector3d & point) { targetVector(point - frame_->position().translation()); }));
}

void VectorOrientationTask::load(mc_solver::QPSolver & solver, const mc_rtc::Configuration & config)
{
  TrajectoryBase::load(solver, config);
  if(config.has("targetVector"))
  {
    targetVector(config("targetVector"));
  }
  if(config.has("relativeVector"))
  {
    Eigen::Matrix3d E_0_r = frame_->robot().posW().rotation().transpose();
    Eigen::Vector3d v = config("relativeVector");
    targetVector(E_0_r * v);
  }
}

} // namespace mc_tasks

namespace
{

static auto registered = mc_tasks::MetaTaskLoader::register_load_function(
    "vectorOrientation",
    [](mc_solver::QPSolver & solver, const mc_rtc::Configuration & config) {
      auto t = [&]() {
        if(config.has("body"))
        {
          return std::make_shared<mc_tasks::VectorOrientationTask>(
              config("body"), config("bodyVector"), solver.robots(),
              robotIndexFromConfig(config, solver.robots(), "vectorOrientation"));
        }
        else
        {
          const auto & robot = robotFromConfig(config, solver.robots(), "vectorOrientation");
          return std::make_shared<mc_tasks::VectorOrientationTask>(robot.frame(config("frame")), config("frameVector"));
        }
      }();
      t->load(solver, config);
      return t;
    });
}

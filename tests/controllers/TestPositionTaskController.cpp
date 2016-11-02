#ifdef BOOST_TEST_MAIN
#undef BOOST_TEST_MAIN
#endif
#include <boost/test/unit_test.hpp>

#include <mc_control/mc_controller.h>
#include <mc_tasks/CoMTask.h>
#include <mc_tasks/PositionTask.h>
#include <mc_control/api.h>

#include <mc_rtc/logging.h>

namespace mc_control
{

struct MC_CONTROL_DLLAPI TestPositionTaskController : public MCController
{
public:
  TestPositionTaskController(std::shared_ptr<mc_rbdyn::RobotModule> rm, double dt)
  : MCController(rm, dt)
  {
    // Check that the default constructor loads the robot + ground environment
    BOOST_CHECK_EQUAL(robots().robots().size(), 2);
    // Check that HRP2-DRC was loaded
    BOOST_CHECK_EQUAL(robot().name(), "hrp2_drc");
    solver().addConstraintSet(contactConstraint);
    solver().addConstraintSet(dynamicsConstraint);
    postureTask->stiffness(1);
    postureTask->weight(1);
    solver().addTask(postureTask.get());
    solver().setContacts({
      mc_rbdyn::Contact(robots(), "LFullSole", "AllGround"),
      mc_rbdyn::Contact(robots(), "RFullSole", "AllGround")
    });

    /* Create and add the position task with the default stiffness/weight */
    posTask = std::make_shared<mc_tasks::PositionTask>("RARM_LINK6", robots(), 0);
    solver().addTask(posTask);

    comTask = std::make_shared<mc_tasks::CoMTask>(robots(), 0);
    solver().addTask(comTask);

    LOG_SUCCESS("Created TestPositionTaskController")
  }

  virtual bool run() override
  {
    bool ret = MCController::run();
    BOOST_CHECK(ret);
    nrIter++;
    if(nrIter == 1999)
    {
      /* Check that the task is "finished" */
      BOOST_CHECK_SMALL(posTask->eval().norm(), 1e-2);
      BOOST_CHECK_SMALL(posTask->speed().norm(), 1e-3);

      /* Apply dimWeight and give a "crazy" position target */
      posTask->dimWeight(Eigen::Vector3d(1., 1., 0.));
      posTask->position(posTask->position() + Eigen::Vector3d(0.1, 0, 100.0));
    }
    if(nrIter == 2999)
    {
      /* At this point the task error on z should be significant */
      BOOST_CHECK(fabs(posTask->eval().z()) > 90);
      /* But the task speed on z should be small */
      BOOST_CHECK_SMALL(fabs(posTask->speed().z()), 1e-3);

      /* Reset the task and ask to raise the hand by 15 cm using only the
       * right arm joints */
      posTask->reset();
      posTask->dimWeight(Eigen::Vector3d(1., 1., 1.));
      posTask->selectActiveJoints(solver(), active_joints);
      posTask->position(posTask->position() + Eigen::Vector3d(0, 0, 0.15));
    }
    if(nrIter == 3999)
    {
      /* Check that the task is "finished" */
      BOOST_CHECK_SMALL(posTask->eval().norm(), 1e-2);
      BOOST_CHECK_SMALL(posTask->speed().norm(), 1e-2);

      /* Now move the hand down again, forbid RAJ3 movement in all tasks */
      posTask->reset();
      posTask->selectUnactiveJoints(solver(), {"RARM_JOINT3"});
      orig_raj3 = robot().mbc().q[robot().jointIndexByName("RARM_JOINT3")][0];
      posTask->position(posTask->position() + Eigen::Vector3d(0, 0, -0.15));

      /* Also reset the joint target in posture task */
      auto p = postureTask->posture();
      p[robot().jointIndexByName("RARM_JOINT3")][0] = orig_raj3;
      postureTask->posture(p);
    }
    if(nrIter == 4999)
    {
      /* Check that the task is "finished" */
      BOOST_CHECK_SMALL(posTask->eval().norm(), 1e-2);
      BOOST_CHECK_SMALL(posTask->speed().norm(), 1e-2);

      /* And that RARM_JOINT3 didn't move. Note that the error is not so
       * small because of other tasks' interaction */
      double current_raj3 = robot().mbc().q[robot().jointIndexByName("RARM_JOINT3")][0];
      BOOST_CHECK_SMALL(fabs(orig_raj3 - current_raj3), 1e-2);
    }
    return ret;
  }

  virtual void reset(const ControllerResetData & reset_data) override
  {
    MCController::reset(reset_data);
    /* Reset the task to the current end-effector position */
    posTask->reset();
    comTask->reset();
    /* Move the end-effector 10cm forward, 10 cm to the right and 10 cm upward */
    posTask->position(posTask->position() + Eigen::Vector3d(0.3, -0.1, 0.2));
  }
private:
  unsigned int nrIter = 0;
  std::shared_ptr<mc_tasks::PositionTask> posTask = nullptr;
  std::shared_ptr<mc_tasks::CoMTask> comTask = nullptr;
  std::vector<std::string> active_joints = [](){
    std::vector<std::string> ret;
    for(unsigned int i = 0; i < 8; ++i)
    {
      std::stringstream ss;
      ss << "RARM_JOINT" << i;
      ret.push_back(ss.str());
    }
    return ret;
  }();
  double orig_raj3 = 0;
};

}

SIMPLE_CONTROLLER_CONSTRUCTOR("TestPositionTaskController", mc_control::TestPositionTaskController)
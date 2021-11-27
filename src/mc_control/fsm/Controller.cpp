/*
 * Copyright 2015-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_control/fsm/Controller.h>

#include <mc_rbdyn/RobotLoader.h>
#include <mc_rbdyn/configuration_io.h>

#include <mc_solver/ConstraintSetLoader.h>

#include <mc_tasks/MetaTaskLoader.h>

#include <mc_rtc/gui/Button.h>
#include <mc_rtc/gui/Form.h>
#include <mc_rtc/gui/Label.h>
#include <mc_rtc/io_utils.h>

namespace mc_rtc
{

mc_control::fsm::Contact ConfigurationLoader<mc_control::fsm::Contact>::load(const mc_rtc::Configuration & config)
{
  return mc_control::fsm::Contact(config("r1"), config("r2"), config("r1Surface"), config("r2Surface"),
                                  config("friction", mc_rbdyn::Contact::defaultFriction),
                                  config("dof", Eigen::Vector6d::Ones().eval()));
}

} // namespace mc_rtc

namespace mc_control
{

namespace fsm
{

namespace
{
/** Always pick a steady clock */
using clock = typename std::conditional<std::chrono::high_resolution_clock::is_steady,
                                        std::chrono::high_resolution_clock,
                                        std::chrono::steady_clock>::type;
} // namespace

#ifdef MC_RTC_BUILD_STATIC
std::unique_ptr<StateFactory> Controller::factory_ptr_;
#endif

Contact Contact::from_mc_rbdyn(const Controller & ctl, const mc_rbdyn::Contact & contact)
{
  return {ctl.robots().robot(contact.r1Index()).name(), ctl.robots().robot(contact.r2Index()).name(),
          contact.r1Surface()->name(), contact.r2Surface()->name(), contact.friction()};
}

Controller::Controller(std::shared_ptr<mc_rbdyn::RobotModule> rm, double dt, const mc_rtc::Configuration & config)
: MCController(std::vector<mc_rbdyn::RobotModulePtr>{rm}, dt, config),
#ifndef MC_RTC_BUILD_STATIC
  factory_(config("StatesLibraries", std::vector<std::string>{}),
           config("StatesFiles", std::vector<std::string>{}),
           config("VerboseStateFactory", false))
#else
  factory_(factory())
#endif
{
#ifdef MC_RTC_BUILD_STATIC
  factory_.load_files(config("StatesFiles", std::vector<std::string>{}));
  factory_.set_verbosity(config("VerboseStateFactory", false));
#endif
  idle_keep_state_ = config("IdleKeepState", false);
  /** Load additional robots from the configuration */
  {
    auto config_robots = config("robots", std::map<std::string, mc_rtc::Configuration>{});
    for(const auto & cr : config_robots)
    {
      const auto & name = cr.first;
      std::string module = cr.second("module");
      auto params = cr.second("params", std::vector<std::string>{});
      mc_rbdyn::RobotModulePtr rm = nullptr;
      if(params.size() == 0)
      {
        rm = mc_rbdyn::RobotLoader::get_robot_module(module);
      }
      else if(params.size() == 1)
      {
        rm = mc_rbdyn::RobotLoader::get_robot_module(module, params.at(0));
      }
      else if(params.size() == 2)
      {
        rm = mc_rbdyn::RobotLoader::get_robot_module(module, params.at(0), params.at(1));
      }
      else
      {
        mc_rtc::log::error_and_throw<std::runtime_error>(
            "FSM controller only handles robot modules that require two parameters at most");
      }
      if(!rm)
      {
        mc_rtc::log::error_and_throw<std::runtime_error>("Failed to load {} as specified in configuration", name);
      }
      loadRobot(rm, name);
      auto & r = robots().robot(name);
      auto & realRobot = realRobots().robot(name);
      r.posW(cr.second("init_pos", sva::PTransformd::Identity()));
      realRobot.posW(r.posW());
    }
    mc_rtc::log::info("Robots loaded in FSM controller:");
    for(const auto & r : robots())
    {
      mc_rtc::log::info("- {}", r.name());
    }
  }
  /** Load global constraints (robots' kinematics/dynamics constraints and contact constraint */
  {
    auto config_constraints = config("constraints", std::vector<mc_rtc::Configuration>{});
    for(const auto & cc : config_constraints)
    {
      constraints_.emplace_back(mc_solver::ConstraintSetLoader::load(solver(), cc));
      if(static_cast<std::string>(cc("type")) == "contact")
      {
        contact_constraint_ = std::static_pointer_cast<mc_solver::ContactConstraint>(constraints_.back());
      }
      /*FIXME Add a name mechanism in ConstraintSet to get information here */
      solver().addConstraintSet(*constraints_.back());
    }
    if(!contact_constraint_)
    {
      mc_rtc::log::warning("No contact constraint loaded from the FSM configuration");
    }
  }
  /** Load collision managers */
  {
    auto config_collisions = config("collisions", std::vector<mc_rtc::Configuration>{});
    for(auto & config_cc : config_collisions)
    {
      if(!config_cc.has("type"))
      {
        config_cc.add("type", "collision");
      }
      auto cc = mc_solver::ConstraintSetLoader::load<mc_solver::CollisionsConstraint>(solver(), config_cc);
      auto & r1 = robots().robot(cc->r1Index);
      auto & r2 = robots().robot(cc->r2Index);
      collision_constraints_[{r1.name(), r2.name()}] = cc;
      solver().addConstraintSet(*cc);
    }
  }
  /** Create posture task for actuated robots */
  for(auto & robot : robots())
  {
    if(robot.mb().nrDof() - robot.mb().joint(0).dof() > 0)
    {
      double stiffness = 1.0;
      double weight = 10.0;
      if(config.has(robot.name()))
      {
        auto robot_config = config(robot.name());
        if(robot_config.has("posture"))
        {
          robot_config("posture")("stiffness", stiffness);
          robot_config("posture")("weight", weight);
        }
      }
      auto t = std::make_shared<mc_tasks::PostureTask>(solver(), robot.robotIndex(), stiffness, weight);
      t->name("FSM_" + t->name());
      posture_tasks_[robot.name()] = t;
      solver().addTask(t);
    }
    if(robot.mb().joint(0).type() == rbd::Joint::Free)
    {
      double stiffness = 2.0;
      double weight = 100.0;
      if(config.has(robot.name()))
      {
        auto robot_config = config(robot.name());
        if(robot_config.has("ff"))
        {
          robot_config("ff")("stiffness", stiffness);
          robot_config("ff")("weight", weight);
        }
      }
      auto t = std::make_shared<mc_tasks::EndEffectorTask>(robot.mb().body(0).name(), solver().robots(),
                                                           robot.robotIndex(), stiffness, weight);
      t->name("FSM_" + t->name());
      ff_tasks_[robot.name()] = t;
    }
  }
  /** Create contacts */
  if(config.has("contacts"))
  {
    contacts_ = config("contacts");
  }
  contacts_changed_ = true;
  /** Load more states if they are provided in the configuration */
  if(config.has("states"))
  {
    factory_.load(config("states"));
  }
  /** Setup executor */
  executor_.init(*this, config_);
}

Controller::~Controller()
{
  executor_.teardown(*this);
  datastore().clear();
}

bool Controller::run()
{
  return run(mc_solver::FeedbackType::None);
}

bool Controller::run(mc_solver::FeedbackType fType)
{
  auto startUpdateContacts = clock::now();
  updateContacts();
  updateContacts_dt_ = clock::now() - startUpdateContacts;
  executor_.run(*this, idle_keep_state_);
  if(!executor_.running())
  {
    if(running_)
    {
      running_ = false;
      startIdleState();
    }
  }
  else
  {
    if(!running_)
    {
      running_ = true;
      teardownIdleState();
    }
  }
  return MCController::run(fType);
}

void Controller::reset(const ControllerResetData & data)
{
  MCController::reset(data);
  if(config().has("init_pos"))
  {
    robot().posW(config()("init_pos"));
    realRobot().posW(robot().posW());
  }
  auto startUpdateContacts = clock::now();
  updateContacts();
  updateContacts_dt_ = clock::now() - startUpdateContacts;

  /** GUI information */
  if(gui_)
  {
    auto all_states = factory_.states();
    std::sort(all_states.begin(), all_states.end());
    gui_->data().add("states", all_states);
    gui_->removeElement({"FSM"}, "Contacts");
    gui_->addElement({"FSM"}, mc_rtc::gui::Label("Contacts", [this]() {
                       std::string ret;
                       for(const auto & c : contacts_)
                       {
                         std::stringstream ss;
                         ss << c.r1Surface << "/" << c.r2Surface << " | " << c.dof.transpose() << " | " << c.friction
                            << "\n";
                         ret += ss.str();
                       }
                       if(ret.size())
                       {
                         ret.pop_back();
                       }
                       return ret;
                     }));
    gui_->removeElement({"Contacts", "Add"}, "Add contact");
    gui_->addElement({"Contacts", "Add"},
                     mc_rtc::gui::Form(
                         "Add contact",
                         [this](const mc_rtc::Configuration & data) {
                           std::string r0 = data("R0");
                           std::string r1 = data("R1");
                           std::string r0Surface = data("R0 surface");
                           std::string r1Surface = data("R1 surface");
                           double friction = data("Friction", mc_rbdyn::Contact::defaultFriction);
                           Eigen::Vector6d dof = data("dof", Eigen::Vector6d::Ones().eval());
                           addContact({r0, r1, r0Surface, r1Surface, friction, dof});
                         },
                         mc_rtc::gui::FormDataComboInput{"R0", true, {"robots"}},
                         mc_rtc::gui::FormDataComboInput{"R0 surface", true, {"surfaces", "$R0"}},
                         mc_rtc::gui::FormDataComboInput{"R1", true, {"robots"}},
                         mc_rtc::gui::FormDataComboInput{"R1 surface", true, {"surfaces", "$R1"}},
                         mc_rtc::gui::FormNumberInput("Friction", false, mc_rbdyn::Contact::defaultFriction),
                         mc_rtc::gui::FormArrayInput<Eigen::Vector6d>("dof", false, Eigen::Vector6d::Ones())));
  }
  logger().addLogEntry("perf_FSM_UpdateContacts", [this]() { return updateContacts_dt_.count(); });
  startIdleState();
}

void Controller::resetPostures()
{
  for(auto & pt : posture_tasks_)
  {
    pt.second->reset();
  }
}

void Controller::startIdleState()
{
  resetPostures();
  // Save posture weights
  saved_posture_weights_.clear();
  for(auto & pt : posture_tasks_)
  {
    saved_posture_weights_[pt.first] = pt.second->weight();
    // Set high weight to prevent the robot from changing configuration
    pt.second->weight(10000);
  }
  for(auto & fft : ff_tasks_)
  {
    fft.second->reset();
    solver().addTask(fft.second);
  }
}

void Controller::teardownIdleState()
{
  // Reset default posture weights
  for(auto & pt : posture_tasks_)
  {
    pt.second->weight(saved_posture_weights_[pt.first]);
  }

  for(auto & fft : ff_tasks_)
  {
    solver().removeTask(fft.second);
  }
}

void Controller::updateContacts()
{
  if(contacts_changed_ && contact_constraint_)
  {
    std::vector<mc_rbdyn::Contact> contacts;
    contact_constraint_->contactConstr->resetDofContacts();

    auto ensureValidContact = [this](const std::string & robotName, const std::string & surfaceName) {
      if(!hasRobot(robotName))
      {
        const auto availableRobots =
            mc_rtc::io::to_string(robots(), [](const mc_rbdyn::Robot & r) { return r.name(); });
        mc_rtc::log::error_and_throw<std::runtime_error>(
            "Failed to add contact: no robot named {} (available robots: {})", robotName, availableRobots);
      }
      if(!robot(robotName).hasSurface(surfaceName))
      {
        mc_rtc::log::error_and_throw<std::runtime_error>("Failed to add contact: no surface named {} in robot {}",
                                                         surfaceName, robotName);
      }
    };
    for(const auto & c : contacts_)
    {
      ensureValidContact(c.r1, c.r1Surface);
      ensureValidContact(c.r2, c.r2Surface);
      auto r1Index = robot(c.r1).robotIndex();
      auto r2Index = robot(c.r2).robotIndex();
      contacts.emplace_back(robots(), r1Index, r2Index, c.r1Surface, c.r2Surface, c.friction);
      auto cId = contacts.back().contactId(robots());
      contact_constraint_->contactConstr->addDofContact(cId, c.dof.asDiagonal());
    }
    solver().setContacts(contacts);
    if(gui_)
    {
      gui_->removeCategory({"Contacts", "Remove"});
      for(const auto & c : contacts_)
      {
        std::string bName = c.r1 + "::" + c.r1Surface + " & " + c.r2 + "::" + c.r2Surface;
        gui_->addElement({"Contacts", "Remove"}, mc_rtc::gui::Button(bName, [this, &c]() { removeContact(c); }));
      }
    }
    contact_constraint_->contactConstr->updateDofContacts();
  }
  contacts_changed_ = false;
}

void Controller::addCollisions(const std::string & r1,
                               const std::string & r2,
                               const std::vector<mc_rbdyn::Collision> & collisions)
{
  if(r1 != r2 && collision_constraints_.count({r2, r1}))
  {
    std::vector<mc_rbdyn::Collision> swapped;
    swapped.reserve(collisions.size());
    for(const auto & c : collisions)
    {
      swapped.push_back({c.body2, c.body1, c.iDist, c.sDist, c.damping});
    }
    addCollisions(r2, r1, swapped);
    return;
  }
  if(!collision_constraints_.count({r1, r2}))
  {
    if(!hasRobot(r1) || !hasRobot(r2))
    {
      mc_rtc::log::error("Try to add collision for robot {} and {} which are not involved in this FSM", r1, r2);
      return;
    }
    auto r1Index = robot(r1).robotIndex();
    auto r2Index = robot(r2).robotIndex();
    collision_constraints_[{r1, r2}] =
        std::make_shared<mc_solver::CollisionsConstraint>(robots(), r1Index, r2Index, solver().dt());
    solver().addConstraintSet(*collision_constraints_[{r1, r2}]);
  }
  auto & cc = collision_constraints_[{r1, r2}];
  mc_rtc::log::info("[FSM] Add collisions {}/{}", r1, r2);
  for(const auto & c : collisions)
  {
    mc_rtc::log::info("[FSM] - {}::{}/{}::{}", r1, c.body1, r2, c.body2);
  }
  cc->addCollisions(solver(), collisions);
}

void Controller::removeCollisions(const std::string & r1,
                                  const std::string & r2,
                                  const std::vector<mc_rbdyn::Collision> & collisions)
{
  if(!collision_constraints_.count({r1, r2}))
  {
    return;
  }
  auto & cc = collision_constraints_[{r1, r2}];
  mc_rtc::log::info("[FSM] Remove collisions {}/{}", r1, r2);
  for(const auto & c : collisions)
  {
    mc_rtc::log::info("[FSM] - {}::{}/{}::{}", r1, c.body1, r2, c.body2);
  }
  cc->removeCollisions(solver(), collisions);
}

void Controller::removeCollisions(const std::string & r1, const std::string & r2)
{
  if(!collision_constraints_.count({r1, r2}))
  {
    return;
  }
  auto & cc = collision_constraints_[{r1, r2}];
  mc_rtc::log::info("[FSM] Remove all collisions {}/{}", r1, r2);
  cc->reset();
}

bool Controller::hasRobot(const std::string & robot) const
{
  return robots().hasRobot(robot);
}

mc_rbdyn::Robot & Controller::robot(const std::string & name)
{
  return robots().robot(name);
}

std::shared_ptr<mc_tasks::PostureTask> Controller::getPostureTask(const std::string & robot)
{
  if(posture_tasks_.count(robot))
  {
    return posture_tasks_.at(robot);
  }
  return nullptr;
}

void Controller::addContact(const Contact & c)
{
  bool inserted;
  ContactSet::iterator it;
  std::tie(it, inserted) = contacts_.insert(c);
  contacts_changed_ |= inserted;
  if(!inserted)
  {
    if(it->dof != c.dof)
    {
      mc_rtc::log::info("[FSM] Changed contact DoF {}::{}/{}::{} to {}", c.r1, c.r1Surface, c.r2, c.r2Surface,
                        c.dof.transpose());
      it->dof = c.dof;
      contacts_changed_ = true;
    }
    if(it->friction != c.friction)
    {
      mc_rtc::log::info("[FSM] Changed contact friction {}::{}/{}::{} to {}", c.r1, c.r1Surface, c.r2, c.r2Surface,
                        c.friction);
      it->friction = c.friction;
      contacts_changed_ = true;
    }
  }
  else
  {
    mc_rtc::log::info("[FSM] Add contact {}::{}/{}::{} (DoF: {})", c.r1, c.r1Surface, c.r2, c.r2Surface,
                      c.dof.transpose());
  }
}

void Controller::removeContact(const Contact & c)
{
  contacts_changed_ |= static_cast<bool>(contacts_.erase(c));
  if(contacts_changed_)
  {
    mc_rtc::log::info("[FSM] Remove contact {}::{}/{}::{}", c.r1, c.r1Surface, c.r2, c.r2Surface);
  }
}

const ContactSet & Controller::contacts() const
{
  return contacts_;
}

bool Controller::hasContact(const Contact & c) const
{
  for(const auto & co : contacts_)
  {
    if(co == c)
    {
      return true;
    }
  }
  return false;
}

bool Controller::resume(const std::string & state)
{
  if(!factory_.hasState(state))
  {
    mc_rtc::log::error("Cannot play unloaded state: {}", state);
    return false;
  }
  return executor_.resume(state);
}

} // namespace fsm

} // namespace mc_control

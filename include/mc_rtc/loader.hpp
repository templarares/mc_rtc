#pragma once

#include <mc_rtc/loader.h>
#include <mc_rtc/loader_sandbox.h>

#include <sstream>

namespace mc_rtc
{

template<typename T>
ObjectLoader<T>::ObjectDeleter::ObjectDeleter(void * sym)
{
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wpedantic"
  delete_fn_ = (void(*)(T*))(sym);
  #pragma GCC diagnostic pop
}

template<typename T>
void ObjectLoader<T>::ObjectDeleter::operator()(T * ptr)
{
  delete_fn_(ptr);
}

template<typename T>
ObjectLoader<T>::ObjectLoader(const std::string & class_name, const std::vector<std::string> & paths, bool enable_sandbox, bool verbose, Loader::callback_t cb)
: class_name(class_name),
  enable_sandbox(enable_sandbox),
  verbose(verbose)
{
  Loader::init();
  load_libraries(paths, cb);
}

template<typename T>
ObjectLoader<T>::~ObjectLoader()
{
  for(const auto & h : handles_)
  {
    lt_dlclose(h.second);
  }
  Loader::close();
}

template<typename T>
bool ObjectLoader<T>::has_object(const std::string & object) const
{
  return handles_.count(object) != 0;
}

template<typename T>
bool ObjectLoader<T>::has_symbol(const std::string & object, const std::string & symbol) const
{
  return has_object(object) && lt_dlsym(handles_.at(object), symbol.c_str()) != nullptr;
}

template<typename T>
std::vector<std::string> ObjectLoader<T>::objects() const
{
  std::vector<std::string> res;
  for(const auto & h : handles_)
  {
    res.push_back(h.first);
  }
  return res;
}

template<typename T>
void ObjectLoader<T>::load_libraries(const std::vector<std::string> & paths,
                                     Loader::callback_t cb)
{
  Loader::load_libraries(class_name, paths, handles_, verbose, cb);
  for(const auto & h : handles_)
  {
    if(deleters_.count(h.first) == 0)
    {
      void * sym = lt_dlsym(h.second, "destroy");
      if(sym == nullptr)
      {
        LOG_ERROR_AND_THROW(LoaderException, "Symbol destroy not found in " << lt_dlgetinfo(h.second)->filename << std::endl << lt_dlerror())
      }
      deleters_[h.first] = ObjectDeleter(sym);
    }
  }
}

template<typename T>
void ObjectLoader<T>::clear()
{
  handles_.clear();
}

template<typename T>
void ObjectLoader<T>::enable_sandboxing(bool enable_sandbox)
{
  this->enable_sandbox = enable_sandbox;
}

template<typename T>
void ObjectLoader<T>::set_verbosity(bool verbose)
{
  this->verbose = verbose;
}

template<typename T>
template<typename... Args>
std::shared_ptr<T> ObjectLoader<T>::create_object(const std::string & name, Args & ... args)
{
  if(!has_object(name))
  {
    LOG_ERROR_AND_THROW(LoaderException, "Requested creation of object named " << name << " which has not been loaded")
  }
  unsigned int args_passed = 1 + sizeof...(Args);
  unsigned int args_required = args_passed;
  void * args_required_sym = lt_dlsym(handles_[name], "create_args_required");
  if(args_required_sym)
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wpedantic"
    auto create_args_required = (unsigned int(*)())(args_required_sym);
    #pragma GCC diagnostic pop
    args_required = create_args_required();
  }
  else
  {
    /* Discard error message */
    lt_dlerror();
  }
  if(args_passed != args_required)
  {
    LOG_ERROR_AND_THROW(LoaderException, args_passed << " arguments passed to create function of " << name << " which expects " << args_required)
  }
  void * sym = lt_dlsym(handles_[name], "create");
  if(sym == nullptr)
  {
    LOG_ERROR_AND_THROW(LoaderException, "Symbol create not found in " << lt_dlgetinfo(handles_[name])->filename << std::endl << lt_dlerror())
  }
  const char * err = lt_dlerror();
  if(err != nullptr)
  {
    LOG_ERROR_AND_THROW(LoaderException, "Failed to resolve create symbol in " << lt_dlgetinfo(handles_[name])->filename << std::endl << err)
  }
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wpedantic"
  std::function<T*(const std::string &, const Args & ...)> create_fn = (T*(*)(const std::string &, const Args & ...))(sym);
  #pragma GCC diagnostic pop
  T * ptr = nullptr;
  if(enable_sandbox)
  {
    ptr = sandbox_function_call(create_fn, name, args...);
  }
  else
  {
    ptr = no_sandbox_function_call(create_fn, name, args...);
  }
  if(ptr == nullptr)
  {
    LOG_ERROR_AND_THROW(LoaderException, "Call to create for object " << name << " failed")
  }
  return std::shared_ptr<T>(ptr, deleters_[name]);
}

} // namespace mc_rtc

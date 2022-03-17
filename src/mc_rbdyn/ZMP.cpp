#include <mc_rbdyn/ZMP.h>
#include <mc_rtc/logging.h>
#include <stdexcept>

namespace mc_rbdyn
{

template<bool ThrowOnSmallForce>
bool zmp(const sva::ForceVecd & netTotalWrench,
                    const Eigen::Vector3d & plane_p,
                    const Eigen::Vector3d & plane_n,
                    Eigen::Vector3d & result,
                    double minimalNetNormalForce)
{
  if(minimalNetNormalForce <= 0)
  {
    mc_rtc::log::error_and_throw("ZMP cannot be computed: the minimalNetNormalForce must be >0 (divide by zero)");
  }

  const Eigen::Vector3d & force = netTotalWrench.force();
  const Eigen::Vector3d & moment_0 = netTotalWrench.couple();
  Eigen::Vector3d moment_p = moment_0 - plane_p.cross(force);
  double floorn_dot_force = plane_n.dot(force);
  // Prevent potential division by zero
  if(floorn_dot_force < minimalNetNormalForce)
  {
    if(ThrowOnSmallForce)
    {
      mc_rtc::log::error_and_throw("ZMP cannot be computed, projected force too small {}", floorn_dot_force);
    }
    mc_rtc::log::error("ZMP cannot be computed, projected force too small {}", floorn_dot_force);
    return false;
  }
  result.noalias() = plane_p + plane_n.cross(moment_p) / floorn_dot_force;
  return true;
}

Eigen::Vector3d zmp(const sva::ForceVecd & netTotalWrench,
                    const Eigen::Vector3d & plane_p,
                    const Eigen::Vector3d & plane_n,
                    double minimalNetNormalForce)
{
  Eigen::Vector3d result;
  zmp<true>(netTotalWrench, plane_p, plane_n, result, minimalNetNormalForce);
  return result;
}

Eigen::Vector3d zmp(const sva::ForceVecd & netWrench, const sva::PTransformd & zmpFrame, double minimalNetNormalForce)
{
  Eigen::Vector3d n = zmpFrame.rotation().row(2);
  Eigen::Vector3d p = zmpFrame.translation();
  return zmp(netWrench, p, n, minimalNetNormalForce);
}

bool maybeZMP(const sva::ForceVecd & netTotalWrench,
              const Eigen::Vector3d & plane_p,
              const Eigen::Vector3d & plane_n,
              Eigen::Vector3d & result,
              double minimalNetNormalForce)
{
  return zmp<false>(netTotalWrench, plane_p, plane_n, result, minimalNetNormalForce);
}

bool maybeZMP(const sva::ForceVecd & netWrench,
              const sva::PTransformd & zmpFrame,
              Eigen::Vector3d & result,
              double minimalNetNormalForce)
{
  Eigen::Vector3d n = zmpFrame.rotation().row(2);
  Eigen::Vector3d p = zmpFrame.translation();
  return maybeZMP(netWrench, p, n, result, minimalNetNormalForce);
}

} // namespace mc_rbdyn

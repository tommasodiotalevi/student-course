#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

using namespace ROOT;
using namespace ROOT::VecOps;

float custom_InvariantMass(RVecF pt, RVecF eta, RVecF phi, RVecF mass){

   const std::size_t size = pt.size();
 
   R__ASSERT(eta.size() == size && phi.size() == size && mass.size() == size);
 
   Float_t x_sum = 0.;
   Float_t y_sum = 0.;
   Float_t z_sum = 0.;
   Float_t e_sum = 0.;
 
   for (std::size_t i = 0u; i < size; ++ i) {
      // Convert to (e, x, y, z) coordinate system and update sums
      const auto x = pt[i] * std::cos(phi[i]);
      x_sum += x;
      const auto y = pt[i] * std::sin(phi[i]);
      y_sum += y;
      const auto z = pt[i] * std::sinh(eta[i]);
      z_sum += z;
      const auto e = std::sqrt(x * x + y * y + z * z + mass[i] * mass[i]);
      e_sum += e;
   }
 
   // Return invariant mass with (+, -, -, -) metric
   return std::sqrt(e_sum * e_sum - x_sum * x_sum - y_sum * y_sum - z_sum * z_sum);

}
#include <iostream>
#include <vector>

#include "ParticleType.hpp"
#include "ResonanceType.hpp"

int main() {
  ParticleType const pt{"pt test", 5.7, -1};

  std::cout << pt.getName() << ' ' << pt.getMass() << ' ' << pt.getCharge()
            << '\n';
  pt.print();

  ResonanceType const rt{"rt test", 1.025, 10, 7};

  std::cout << rt.getName() << ' ' << rt.getMass() << ' ' << rt.getCharge()
            << ' ' << rt.getWidth() << '\n';
  rt.print();

  std::vector<const ParticleType*> ptr_vector{};
  ptr_vector.push_back(&pt);
  ptr_vector.push_back(&rt);
  auto ptr_vector_size = ptr_vector.size();

  for (auto const& p : ptr_vector) {
    p->print();
  }
}

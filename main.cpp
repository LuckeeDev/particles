#include <iostream>
#include <vector>

#include "Particle.hpp"
#include "ParticleType.hpp"
#include "ResonanceType.hpp"

int main() {
  std::cout << "TESTING THE \"ParticleType\" AND \"ResonanceType\" CLASSES"
            << '\n';

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

  for (auto const& p : ptr_vector) {
    p->print();
  }

  std::cout << "\n\n"
            << "TESTING THE \"Particle\" CLASS" << '\n';

  Particle::addParticleType("electron", 5.7, -1);
  Particle::addParticleType("proton", 1.025, 10, 7);

  Particle e{"electron"};
  Particle p{"proton"};
  Particle n{"neutron"};

  e.printData();
  p.printData();
  n.printData();

  Particle::addParticleType("electron", 1.025, 0, 7);
  Particle::addParticleType("neutron", 1.025, 0, 7);
  n.setIndex("neutron");
  n.printData();

  Particle::printParticleTypes();
}

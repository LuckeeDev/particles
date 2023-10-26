#include "Particle.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "ParticleType.hpp"
#include "ResonanceType.hpp"

// init static members

std::vector<ParticleType> Particle::m_particle_types{};

// momentum operators

Momentum Momentum::operator+(Momentum const& momentum) const {
  return {x + momentum.x, y + momentum.y, z + momentum.z};
}

double Momentum::operator*(Momentum const& momentum) const {
  return x * momentum.x + y * momentum.y + z * momentum.z;
}

// constructor

Particle::Particle(std::string const& name, Momentum const& momentum)
    : m_momentum{momentum} {
  m_index = mFindParticleIndex(name);

  if (m_index == std::nullopt) {
    std::cout << "The \"" << name << "\" particle type does not exist!" << '\n';
  }
}

// public methods

void Particle::printData() const {
  if (m_index != std::nullopt) {
    std::cout << "Index: " << m_index.value() << '\n'
              << "Name: " << m_particle_types[m_index.value()].getName() << '\n'
              << "Momentum: (" << m_momentum.x << ", " << m_momentum.y << ", "
              << m_momentum.z << ")\n";
  }
}

// getters

std::optional<int> Particle::getIndex() const { return m_index; }

Momentum Particle::getMomentum() const { return m_momentum; }

double Particle::getEnergy() const {
  if (m_index != std::nullopt) {
    return std::sqrt(std::pow(m_particle_types[m_index.value()].getMass(), 2) +
                     m_momentum * m_momentum);
  } else {
    std::cout << "This particle has no energy because its index is invalid!"
              << '\n';

    return 0;
  }
}

double Particle::getMass() const {
  if (m_index) {
    return m_particle_types[m_index.value()].getMass();
  } else {
    std::cout << "This particle has no mass because its index is invalid!"
              << '\n';

    return 0;
  }
}

double Particle::getInvariantMass(Particle const& p) const {
  auto newMomentum = m_momentum + p.getMomentum();
  return std::sqrt(std::pow(getEnergy() + p.getEnergy(), 2) -
                   newMomentum * newMomentum);
}

// setters

void Particle::setIndex(std::string const& name) {
  m_index = mFindParticleIndex(name);
}

void Particle::setIndex(int index) {
  if (index >= 0 && index < static_cast<int>(m_particle_types.size())) {
    m_index = index;
  } else {
    std::cout << "The index \"" << index
              << "\" does not refer to a particle type!" << '\n';
  }
}

void Particle::setMomentum(double px, double py, double pz) {
  m_momentum = {px, py, pz};
}

void Particle::setMomentum(Momentum const& momentum) { m_momentum = momentum; }

// static methods

int Particle::countParticleTypes() { return m_particle_types.size(); }

void Particle::addParticleType(std::string const& name, double mass, int charge,
                               double width) {
  auto existing_index = mFindParticleIndex(name);

  if (existing_index == std::nullopt) {
    if (width == 0) {
      m_particle_types.push_back(ParticleType{name, mass, charge});
    } else {
      m_particle_types.push_back(ResonanceType{name, mass, charge, width});
    }
  } else {
    std::cout << "The \"" << name << "\" particle type already exists!" << '\n';
  }
}

void Particle::printParticleTypes() {
  for (auto const& p : m_particle_types) {
    p.print();
    std::cout << '\n';
  }
}

// private methods

std::optional<int> Particle::mFindParticleIndex(std::string const& name) {
  auto v_begin = m_particle_types.begin();
  auto v_end = m_particle_types.end();

  auto it = std::find_if(v_begin, v_end, [&name](ParticleType const& pt) {
    return pt.getName() == name;
  });

  if (it == v_end) {
    return std::nullopt;
  }

  return std::distance(m_particle_types.begin(), it);
}
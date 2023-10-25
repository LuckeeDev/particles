#include "Particle.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>

#include "ParticleType.hpp"
#include "ResonanceType.hpp"

// momentum operators

Momentum Momentum::operator+(Momentum const& momentum) const {
  return {x + momentum.x, y + momentum.y, z + momentum.z};
}

double Momentum::operator*(Momentum const& momentum) const {
  return x * momentum.x + y * momentum.y + z * momentum.z;
}

// constructor

Particle::Particle(std::string const& name,
                   Momentum const& momentum = {0, 0, 0})
    : m_momentum{momentum} {
  m_index = mFindParticleIndex(name);
}

// public methods

void Particle::printData() const {
  std::cout << "Index :" << m_index << '\n'
            << "Name: " << m_particle_types[m_index]->getName() << '\n'
            << "Momentum: (" << m_momentum.x << ", " << m_momentum.y << ", "
            << m_momentum.z << ")\n";
}

// getters

int Particle::getIndex() const { return m_index; }

Momentum Particle::getMomentum() const { return m_momentum; }

double Particle::getEnergy() const {
  return std::sqrt(std::pow(m_particle_types[m_index]->getMass(), 2) +
                   m_momentum * m_momentum);
}

double Particle::getMass() const {
  return m_particle_types[m_index]->getMass();
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
  m_index = m_particle_types[index] == nullptr ? -1 : index;
}

void Particle::setMomentum(double px, double py, double pz) {
  m_momentum = {px, py, pz};
}

void Particle::setMomentum(Momentum const& momentum) { m_momentum = momentum; }

// static methods

int Particle::getParticleTypesCount() { return m_particle_types.size(); }

void Particle::addParticleType(std::string const& name, double mass, int charge,
                               double width = 0) {
  auto existing_index = mFindParticleIndex(name);

  if (existing_index != -1) {
    auto first_empty =
        std::find(m_particle_types.begin(), m_particle_types.end(), nullptr);

    ParticleType* new_particle;

    if (width == 0) {
      new_particle = new ParticleType{name, mass, charge};
    } else {
      new_particle = new ResonanceType{name, mass, charge, width};
    }

    *first_empty = new_particle;
  }
}

void Particle::printParticleTypes() {
  for (auto const& p : m_particle_types) {
    p->print();
    std::cout << '\n';
  }
}

// private methods

int Particle::mFindParticleIndex(std::string const& name) {
  auto it =
      std::find_if(m_particle_types.begin(), m_particle_types.end(),
                [&name](ParticleType* pt) { return pt->getName() == name; });

  if (it == m_particle_types.end()) {
    std::cout << "Invalid particle name!" << '\n';
    return -1;
  }

  return std::distance(m_particle_types.begin(), it);
}
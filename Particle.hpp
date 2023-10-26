#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <array>
#include <optional>
#include <string>

class ParticleType;

struct Momentum {
  double x;
  double y;
  double z;

  Momentum operator+(Momentum const&) const;
  double operator*(Momentum const&) const;
};

class Particle {
 public:
  Particle(std::string const&, Momentum const& = {0., 0., 0.});

  void printData() const;

  // setters

  void setIndex(std::string const&);
  void setIndex(int);
  void setMomentum(double, double, double);
  void setMomentum(Momentum const&);

  // getters

  int getIndex() const;
  Momentum getMomentum() const;
  double getEnergy() const;
  double getMass() const;
  double getInvariantMass(Particle const&) const;

  // static methods

  static int countParticleTypes();
  static void addParticleType(std::string const&, double, int, double = 0.);
  static void printParticleTypes();

 private:
  Momentum m_momentum;
  int m_index;

  static std::array<std::optional<ParticleType>, 10> m_particle_types;
  static int mFindParticleIndex(std::string const& name);
};

#endif
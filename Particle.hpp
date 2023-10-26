#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <optional>
#include <string>
#include <vector>

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

  std::optional<int> getIndex() const;
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
  std::optional<int> m_index;

  static std::vector<ParticleType> m_particle_types;
  static std::optional<int> mFindParticleIndex(std::string const& name);
};

#endif
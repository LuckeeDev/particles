#include <string>

#ifndef PARTICLE_TYPE_HPP
#define PARTICLE_TYPE_HPP

class ParticleType {
 public:
  ParticleType(std::string const&, double, int);
  virtual ~ParticleType() = default;

  std::string getName() const;
  double getMass() const;
  int getCharge() const;
  virtual void print() const;

 private:
  std::string m_name;
  double m_mass;
  int m_charge;
};

#endif
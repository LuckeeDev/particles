#ifndef RESONANCE_TYPE_HPP
#define RESONANCE_TYPE_HPP

#include "ParticleType.hpp"

class ResonanceType : public ParticleType {
 public:
  ResonanceType(std::string const&, double, int, double);

  double getWidth() const;
  void print() const override;

 private:
  double const m_width;
};

#endif
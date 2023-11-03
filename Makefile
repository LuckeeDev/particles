root:
	root -l -b -q -e '.L ParticleType.cpp++'
	root -l -b -q -e '.L ResonanceType.cpp++'
	root -l -b -q -e '.L Particle.cpp++'
	root -e 'gROOT->LoadMacro("generate.cpp")'
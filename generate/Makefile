root:
	root -l -b -q -e '.L ParticleType.cpp++'
	root -l -b -q -e '.L ResonanceType.cpp++'
	root -l -b -q -e '.L Particle.cpp++'
	root -e 'gROOT->LoadMacro("generate.cpp")'

test:
	g++ ParticleType.cpp ResonanceType.cpp Particle.cpp test_main.cpp `root-config --glibs --cflags --libs` -o particles_test.out
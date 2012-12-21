lattice_qcd
===========

references:

	GPU Implementation of the Feynman Path-Integral Method in Quantum Mechanics
            publications.lib.chalmers.se/records/fulltext/144105.pdf

	Lattice QCD for novices 
	    http://arxiv.org/abs/hep-lat/0506036

2. lattice QCD for people who want results http://arxiv.org/pdf/hep-lat/0509046.pdf
3. QCD on the lattice http://www.springerlink.com/content/978-3-642-01849-7/contents/
4. http://scholar.google.com/scholar?q=lattice+qcd&btnG=&hl=en&as_sdt=0%2C22
5. rev.mod.physics http://scholar.google.com/scholar?q=lattice+qcd&btnG=&hl=en&as_sdt=0%2C22
5. rev.mod.physics http://goo.gl/dODgR
HTTP://RMP.APS.ORG/ABSTRACT/RMP/V84/I2/P449_1
HTTP://RMP.APS.ORG/ABSTRACT/RMP/V67/I4/P893_1
6. LATTICE QCD ON PCS? HTTP://WWW.SCIENCEDIRECT.COM/SCIENCE/ARTICLE/PII/S0920563201016395
7. what is renormalization? http://arxiv.org/pdf/hep-ph/0506330v1.pdf



TODO:
	* FIND OUT PARAMETERS RELATION:
       	    It seems that if you increase N, you also have to increase nrepeat so the paths are not correlated and it works. So, either (N=10,nrepeat=10) or (N=50,nrepeat=100) works well. Using a=10 for now.
       	    A nice graph about parameters can be found on GPU Implem p. 75 (it relates a to N)

	* Implement error function (bin,bootstraps)

	Normalization
	    why 1/(2pi*a)**(N/2) doesnt work?

	Calculate acceptance rate (so we can tune epsilon)


todorht:

1. fermion lagrangian? http://www-zeuthen.desy.de/students/2006/doc/nube-marinkovic.pdf
2. add interaction terms
3. implement feynman-kleinert formula
4. QCD phase diagram
5. wtf is Brodsky-Lepage-Mackenzie?
6. stalk lepage's publication


thoughts:
1. heisenberg picture correlation functions
2. 



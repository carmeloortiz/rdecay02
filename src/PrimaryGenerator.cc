/// \file PrimaryGenerator.cc
/// \brief Implementation of the PrimaryGenerator1 class

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <random>
#include <math.h>

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include "PrimaryGenerator.hh"
#include "CEvNSGlow.hh"


PrimaryGenerator::PrimaryGenerator()
  : G4VPrimaryGenerator(), fFSNeutrino(false)
{ }

PrimaryGenerator::~PrimaryGenerator()
{ }

//..............................................................................

void PrimaryGenerator::GeneratePrimaryVertexOpt(G4Event* event, std::vector<double> &xyzb)
{
  const G4int n_particle = 1250; // from our DkMatter paper, via SCENE. 1250 photons per 100 keV n.r.

  const G4double x = xyzb.at(0)*(G4UniformRand()-0.5) ;  
  const G4double y = xyzb.at(1)*(G4UniformRand()-0.5) ;  
  const G4double z = xyzb.at(2)*(G4UniformRand()-0.5) ; 
  G4ThreeVector positionA(x,y,z);

  for (int ii =0; ii< n_particle; ii++) 
    {
      G4double alpha = pi*2.*(G4UniformRand()-0.5);     //alpha uniform in (-pi, pi)
      G4double beta = pi*G4UniformRand();     //alpha uniform in (0, pi)
      G4double ux = std::cos(alpha)*std::sin(beta);
      G4double uy = std::sin(alpha)*std::sin(beta);
      G4double uz = std::cos(beta);

      G4double timeA = 0*s;
      G4PrimaryVertex* vertexA = new G4PrimaryVertex(positionA, timeA);
  
      //primary particle at vertexA
      G4ParticleDefinition* particleDefinition
	= G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");
      G4PrimaryParticle* particle1 = new G4PrimaryParticle(particleDefinition);
      particle1->SetMomentumDirection(G4ThreeVector(ux,uy,uz));    
      particle1->SetKineticEnergy(9.686 * eV);

      vertexA->SetPrimary(particle1);
      event->AddPrimaryVertex(vertexA);
    }

  SetParticlePosition(positionA);
}

//..............................................................................

void PrimaryGenerator::GeneratePrimaryVertexMarley(G4Event* event, std::vector<double> &xyzb, marley::Event marlev) 
{
  // Loop over each of the final particles in the MARLEY event
  const G4double x = xyzb.at(0) *2*(G4UniformRand()-0.5) ;
  const G4double y = xyzb.at(1) *2*(G4UniformRand()-0.5) ;
  const G4double z = xyzb.at(2) *2*(G4UniformRand()-0.5) ;

  const G4int n_particle = 1250; // from our DkMatter paper, via SCENE. 1250 photons per 100 keV n.r.                                                           

  G4ThreeVector positionA(x,y,z);
  G4double timeA = 0*s;
  G4PrimaryVertex* vertex = new G4PrimaryVertex(positionA, timeA);
  fFSNeutrino = false;
      
  for ( const auto& fp : marlev.get_final_particles() ) {

    // Convert each one from a marley::Particle into a G4PrimaryParticle.
    // Do this by first setting the PDG code and the 4-momentum components.
    // If this is an Ar nuclear recoil, let's place opticalphotons on stack, instead.

    if( fp->pdg_code() == 12)
      fFSNeutrino = true;

    if (fp->pdg_code() == 1000180400) {
      int Nphotons_stack = fp->kinetic_energy() * 1000. * n_particle/100.;
      G4ParticleDefinition* particleDefinition
	= G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");

      for (int phot=0; phot<Nphotons_stack; phot++)
	{
	  G4double alpha = pi*2.*(G4UniformRand()-0.5);     //alpha uniform in (-pi, pi)
	  G4double beta = pi*G4UniformRand();     //alpha uniform in (0, pi)
	  G4double ux = std::cos(alpha)*std::sin(beta);
	  G4double uy = std::sin(alpha)*std::sin(beta);
	  G4double uz = std::cos(beta);

	  G4PrimaryParticle* particle = new G4PrimaryParticle(particleDefinition);
	  particle->SetMomentumDirection(G4ThreeVector(ux,uy,uz));
	  particle->SetKineticEnergy(9.686 * eV); // 128nm
	  vertex->SetPrimary(particle);
	}
      if (Nphotons_stack)
	std::cout << "GeneratePrimaryVertexMarley: Put " << int(Nphotons_stack) << " 9.7 eV photons on stack in lieu of Ar nucleus of KE " << fp->kinetic_energy()*1000. << " keV" << std::endl;

    }
    else {
      G4PrimaryParticle* particle = new G4PrimaryParticle( fp->pdg_code(),
							   fp->px(), fp->py(), fp->pz(), fp->total_energy() );

      // Also set the charge of the G4PrimaryParticle appropriately
      particle->SetCharge( fp->charge() );

      // Add the fully-initialized G4PrimaryParticle to the primary vertex
      vertex->SetPrimary( particle );
    }

  } // end of loop on particles

  // The primary vertex has been fully populated with all final-state particles
  // from the MARLEY event. Add it to the G4Event object so that Geant4 can
  // begin tracking the particles through the simulated geometry.

  SetParticlePosition(positionA);
  event->AddPrimaryVertex( vertex );

}

//..............................................................................
// (EDIT): NEW FUNCTION GeneratePrimaryVertexCEvNSGlow
// Repeatedly reads events.dat, and generates primaries according to the contents of the file

void PrimaryGenerator::GeneratePrimaryVertexCEvNSGlow(G4Event* event, std::vector<double> &xyzb) 
{
	double time;

	//opens index.dat (file which saves the index of the current event)
	std::fstream index_file;
	index_file.open("CEvNSGlow/index.dat", std::fstream::in);
	
	//opens events.dat
	std::fstream events_file;
	events_file.open("CEvNSGlow/events.dat", std::fstream::in);

	//opens hits.dat
	std::fstream hits_file;
	hits_file.open("CEvNSGlow/hits.dat", std::fstream::app);

	//reads the index, and chooses the line corresponding to the index
	std::string index;
	std::string line;
	std::getline(index_file,index);

	int i = std::stoi(index);
	while(i>=0){
		std::getline(events_file,line);
		i--;
	}
	
	//runs different things according to what event is to be run
	int tag;
	int j;
	j = line.find(' '); 
	tag = std::stoi(line.substr(0, j));

	if(tag==0){
		//std::cout << "Running a CEvNS event!" << std::endl;
		fFSNeutrino = false;
		hits_file << 0 << std::endl;

		//position and time initialization
		const G4double x = xyzb.at(0) *2*(G4UniformRand()-0.5) ;
		const G4double y = xyzb.at(1) *2*(G4UniformRand()-0.5) ;
		const G4double z = xyzb.at(2) *2*(G4UniformRand()-0.5) ;
		double E_rec;
		int i_1;
		int i_2;
		int n_particle = 1250;

		i_1 = line.find(' '); 
		i_2 = line.find(' ',i_1+1);
		time = stod(line.substr(i_1+1, i_2));
	
		G4ThreeVector positionA(x,y,z);
		G4double timeA = time*s;
		G4PrimaryVertex* vertex = new G4PrimaryVertex(positionA, timeA);

		//energy
		i_1 = i_2; 
		i_2 = line.find(' ',i_1+1);
		E_rec = stod(line.substr(i_1+1, 10));
	
		//particle declaration and assignment
		int Nphotons_stack = E_rec* 1000. * n_particle/100.;
	      	G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");

	      	for (int phot=0; phot<Nphotons_stack; phot++){
			G4double alpha = pi*2.*(G4UniformRand()-0.5);     //alpha uniform in (-pi, pi)
			G4double beta = pi*G4UniformRand();     //alpha uniform in (0, pi)				
			G4double ux = std::cos(alpha)*std::sin(beta);
			G4double uy = std::sin(alpha)*std::sin(beta);
			G4double uz = std::cos(beta);

			G4PrimaryParticle* particle = new G4PrimaryParticle(particleDefinition);
			particle->SetMomentumDirection(G4ThreeVector(ux,uy,uz));
			particle->SetKineticEnergy(9.686 * eV); // 128nm
			vertex->SetPrimary(particle);
		}
		if (Nphotons_stack){
			std::cout << "GeneratePrimaryVertexMarley: Put " << int(Nphotons_stack);
			std::cout << " 9.7 eV photons on stack in lieu of Ar nucleus of KE ";
			std::cout << E_rec*1000. << " keV" << std::endl;
	 	}
		
		SetParticlePosition(positionA);
	  	event->AddPrimaryVertex( vertex );
		
		std::cout << "PARTICLE DETECTED " << particleDefinition->GetParticleName() << std::endl;
		
	}
	else if(tag==1){
		
		//std::cout << "Running a CC event!" << std::endl;
		//For CC events we use MARLEY
		hits_file << 1 << std::endl;

		//Collect the beta value and E_ave from line
		double beta;
		double E_ave;
		int i_1;
		int i_2;

		//time
		i_1 = line.find(' '); 
		i_2 = line.find(' ',i_1+1);
		time = stod(line.substr(i_1+1, i_2-i_1));

		//beta
		i_1 = i_2; 
		i_2 = line.find(' ',i_1+1);
		beta = std::stod(line.substr(i_1+1, i_2-i_1));

		//E_ave
		i_1 = i_2; 
		i_2 = line.find(' ',i_1+1);
		E_ave = std::stod(line.substr(i_1+1, 10));

		std::cout << "time: " << time << " s" <<std::endl;
		std::cout << "beta: " << beta <<std::endl;
		std::cout << "E_ave: " << E_ave << " MeV" <<std::endl;

		//Generate Marley Event
		std::string config_file_name = "/home/local1/Downloads/rdecay02-liquid_deception/build/CEvNSGlow/CEvNS.js";
		marley::JSONConfig config(config_file_name);
		marley_generator_= config.create_generator(E_ave, beta); //E_ave, beta value
          	marley::Event marlev = marley_generator_.create_event();

		//Feed the products into Geant4

		//position and time initialization
		const G4double x = xyzb.at(0) *2*(G4UniformRand()-0.5) ;
		const G4double y = xyzb.at(1) *2*(G4UniformRand()-0.5) ;
		const G4double z = xyzb.at(2) *2*(G4UniformRand()-0.5) ;
		const G4int n_particle = 1250;  
		
		G4ThreeVector positionA(x,y,z);
		G4double timeA = time*s;
		G4PrimaryVertex* vertex = new G4PrimaryVertex(positionA, timeA);
		fFSNeutrino = false;

		// Loop over each of the final particles in the MARLEY event    
	  	for ( const auto& fp : marlev.get_final_particles() ) {

		if( fp->pdg_code() == 12)
			fFSNeutrino = true;

	        if (fp->pdg_code() == 1000180400) {
	      		int Nphotons_stack = fp->kinetic_energy() * 1000. * n_particle/100.;
	      		G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");

	      		for (int phot=0; phot<Nphotons_stack; phot++){
				G4double alpha = pi*2.*(G4UniformRand()-0.5);     //alpha uniform in (-pi, pi)
				G4double beta = pi*G4UniformRand();     //alpha uniform in (0, pi)
				G4double ux = std::cos(alpha)*std::sin(beta);
				G4double uy = std::sin(alpha)*std::sin(beta);
				G4double uz = std::cos(beta);

				G4PrimaryParticle* particle = new G4PrimaryParticle(particleDefinition);
				particle->SetMomentumDirection(G4ThreeVector(ux,uy,uz));
				particle->SetKineticEnergy(9.686 * eV); // 128nm
				vertex->SetPrimary(particle);
			}
			if (Nphotons_stack){
				std::cout << "GeneratePrimaryVertexMarley: Put " << int(Nphotons_stack);
				std::cout << " 9.7 eV photons on stack in lieu of Ar nucleus of KE ";
				std::cout << fp->kinetic_energy()*1000. << " keV" << std::endl;
		 	}
	        }
	        else {
		        G4PrimaryParticle* particle = new G4PrimaryParticle( fp->pdg_code(),fp->px(), fp->py(), fp->pz(), fp->total_energy() );
	      		particle->SetCharge( fp->charge() );
	      		vertex->SetPrimary( particle );
	        }

		} //end of for loop on MARLEY event final particles
		SetParticlePosition(positionA);
	  	event->AddPrimaryVertex( vertex );
		
	}
	else if(tag==2){
		
		//std::cout << "Running an Ar39 decay event!" << std::endl;
		//Ar39 decay
		fFSNeutrino = false;
		hits_file << 2 << std::endl;

		//position and time initialization
		const G4double x = xyzb.at(0) *2*(G4UniformRand()-0.5) ;
		const G4double y = xyzb.at(1) *2*(G4UniformRand()-0.5) ;
		const G4double z = xyzb.at(2) *2*(G4UniformRand()-0.5) ;
		int i_1;
		int i_2;

		i_1 = line.find(' '); 
		i_2 = line.find(' ',i_1+1);
		time = stod(line.substr(i_1+1, 10));
	
		G4ThreeVector positionA(x,y,z);//positionA(x,y,z);
		G4double timeA = time*s;
		G4PrimaryVertex* vertex = new G4PrimaryVertex(positionA, 0); //timeA

		//TEMPORARY EDIT
		std::fstream out_events;
		out_events.open("CEvNSGlow/Ar39_position.dat",std::fstream::app);
		out_events << x << " " << y << " " << z << std::endl;
		out_events.close();
		//TEMPORARY EDIT

		//momentum
		
		G4double alpha = pi*2.*(G4UniformRand()-0.5);     //alpha uniform in (-pi, pi)
		G4double beta = pi*G4UniformRand();     //alpha uniform in (0, pi)
		G4double ux = std::cos(alpha)*std::sin(beta);
		G4double uy = std::sin(alpha)*std::sin(beta);
		G4double uz = std::cos(beta);
	
		//particle declaration and assignment
		G4IonTable IonTable;
		//std::cout << IonTable.GetNumberOfElements() << std::endl;
		//if (std::stoi(index)==0){
		//	IonTable.CreateIon(18,42,0);
		//}
		//std::cout << IonTable.GetNumberOfElements() << std::endl;
	        G4ParticleDefinition* particleDefinition = IonTable.GetIon(18,39,0,0);//G4IonTable::GetIonTable()->GetIon(18,39,0,0);
		//particleDefinition->SetPDGStable(false);
		//particleDefinition->SetPDGLifeTime(0*second);
		
		//IonTable.GetIon(18,42,0,0);
		//G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");
		G4PrimaryParticle* particle = new G4PrimaryParticle(particleDefinition);
		particle->SetProperTime(1*ns);
		particle->SetMomentumDirection(G4ThreeVector(ux,uy,uz));
		particle->SetKineticEnergy(0.0*eV);
		particle->SetCharge(0);
		vertex->SetPrimary(particle);
		
		SetParticlePosition(positionA);
	  	event->AddPrimaryVertex( vertex );
		
		std::cout << "PARTICLE DETECTED " << particleDefinition->GetParticleName() << std::endl;
		
		
	}
	index_file.close();
	index_file.open("CEvNSGlow/index.dat", std::fstream::out);

	index_file << std::stoi(index)+1;
	index_file.close();

	hits_file << std::setprecision(10) << time << std::endl;
	hits_file.close();
	events_file.close();
}




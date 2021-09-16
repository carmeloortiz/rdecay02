{
  //seed: , // Random number seed (omit to use time since Unix epoch)

  // Pure 40Ar target
  target: {
    nuclides: [ 1000180400 ],
    atom_fractions: [ 1.0 ],
  },

  // Simulate CC ve scattering and CEvNS on 40Ar
  reactions: ["CEvNS40Ar.react"],//"ve40ArCC_Bhattacharya2009.react", "ES.react","CEvNS40Ar.react"

  //source
  //FOR LATER: TH1 source
  //source: {
  //  type: "th1",
  //  neutrino: "ve",
  //  tfile: "nu-energy-spectrum.root",  // Name of the ROOT file containing
  //                               // the TGraph object

  //  namecycle: "nu-energy-spectrum",        // Name of the TGraph object (used to
  //                               // retrieve it from the ROOT file)
  //},

  //"BETA FIT"
  //
  source: {
      type: "beta-fit",
      neutrino: "ve",
      Emin: 0,           // Minimum neutrino energy (MeV)
      Emax: 60,          // Maximum neutrino energy (MeV)
      Emean:10,         // Mean neutrino energy (MeV)
      beta: 5.0,         // Pinching parameter (dimensionless, default 4.5)
    },

  //source: {
  //  type: "fermi-dirac",
  //  neutrino: "ve",
  //  Emin: 0,           // Minimum neutrino energy (MeV)
  //  Emax: 60,          // Maximum neutrino energy (MeV)
  //  temperature: 3.5,  // Temperature (MeV)
  //  eta: 4             // Pinching parameter (dimensionless, default 0)
  //},

//  source: {
//    neutrino: "ve",        // The source produces electron neutrinos
//    type: "monoenergetic",
//    energy: 20.0,          // MeV
//  },

  // Incident neutrino direction 3-vector
  direction: { x: 0.0, y: 0.0, z: 1.0 },

  // Settings for marley command-line executable
  executable_settings: {

    // The number of events to generate
    events: 1,

    // Event output configuration
    output: [ { file: "events.ascii", format: "ascii", mode: "overwrite" } ],

  },
}

/*
 */

//#include "NangaParbat/preprocessing.h"

#include <yaml-cpp/yaml.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <utility>   // std::pair
#include <stdexcept> // std::runtime_error
#include <sstream>   // std::stringstream

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

// Structure to hold the data for each row
struct DataRow {
    int index;
    double Eb;
    double x, y, z, Q2, pT;
    double xl,xr,zl,zr,Q2l,Q2r,pTl,pTr;
    double value, stat;
    std::string had;
};
// Function to parse a line of CSV data and return a DataRow
DataRow parseLine(const std::string& line) {
    std::stringstream ss(line);
    std::string item;
    DataRow row;

    // Read the index
    std::getline(ss, item, ',');
    row.index = std::stoi(item);

    std::getline(ss, item, ','); 
    row.Eb = std::stod(item);
    // Read x
    std::getline(ss, item, ',');
    if (std::stod(item)==0) row.x=0.00001;
    else row.x = std::stod(item);


    // next 2 columns (xl,xr)
    std::getline(ss, item, ','); // xl
    row.xl=std::stod(item);
    std::getline(ss, item, ','); // xr
    row.xr=std::stod(item);

    // Read y
    std::getline(ss, item, ',');
    row.y = std::stod(item);

    // Read z
    std::getline(ss, item, ',');
    row.z = std::stod(item);

    //  next 2 columns (zl, zr)
    std::getline(ss, item, ','); // zl
    row.zl=std::stod(item);
    std::getline(ss, item, ','); // zr
    row.zr=std::stod(item);

    // Read Q2
    std::getline(ss, item, ',');
    row.Q2 = std::stod(item);

    // next 2 columns (Q2l, Q2r)
    std::getline(ss, item, ','); // Q2l
    row.Q2l=std::stod(item);
    std::getline(ss, item, ','); // Q2r
    row.Q2r=std::stod(item);

    // Read pT
    std::getline(ss, item, ','); 
    row.pT = std::stod(item);

    // next 3 columns (pTl, pTr, obs)
    std::getline(ss, item, ','); // pTl
    row.pTl=std::stod(item);
    std::getline(ss, item, ','); //  pTr
    row.pTr=std::stod(item);
    std::getline(ss, item, ','); // Skip obs

    // Read value
    std::getline(ss, item, ',');
    row.value = std::stod(item);

    // Read stat
    std::getline(ss, item, ',');
    row.stat = std::stod(item);

    // Skip the next 4 columns (syst, target, hadron)
    std::getline(ss, item, ','); // Skip syst
    std::getline(ss, item, ','); // Skip target
    std::getline(ss, item, ','); // hadron
    row.had = item;
    // Return the populated DataRow
    return row;
}


namespace NangaParbat
{

  //_________________________________________________________________________________
  std::string PreprocessSoLID(std::string const& RawDataPath, std::string const& ProcessedDataPath)
  {
    std::cout << "Processing SoLID Pion projection ..." << std::endl;

    // Path to the raw-data folder
    const std::string RawDataFolder = RawDataPath + "/SoLID_Unpol/";


    // Vector of tables to process
    const std::vector<std::string> tables = {"numbers_x1_pip_11_He3.txt","numbers_x1_pim_11_He3.txt","numbers_x2_pip_11_He3.txt","numbers_x2_pim_11_He3.txt","numbers_hQ2_x2_pip_11_He3.txt","numbers_hQ2_x2_pim_11_He3.txt","numbers_x1_pip_11_Q22.5-3_He3.txt","numbers_x1_pim_11_Q22.5-3_He3.txt"};

    // Output folder
    const std::string ofolder = "oSoLID_Unpol";

    // Vector for output names
    std::vector<std::string> filenames;

    // Loop over tables
    for (auto const& tab : tables)
    {
      //naming map
      std::map<std::string, std::pair<double, double>> Q2rangelims;
      std::map<std::string, std::pair<double, double>> xrangelims;
      std::map<std::string, std::pair<double, double>> zrangelims;
      
      //std::cout<<"check names "<<tab<<std::endl;
      // Vs ( = sqrt(2*M*E) for fixed target experiments)
      double Vs;
      if(tab=="numbers_x1_pip_11_He3.txt" or tab=="numbers_x1_pim_11_He3.txt"){
        //std::cout<<" check condition 1"<<std::endl;
        // Initialize naming map for the Q-integration ranges (the first element is for the name of the output data file)
        Q2rangelims = {{"_11_Q2_1_1.50", {1.0, 1.50}}, {"_11_Q2_1.50_2.0", {1.50, 2.0}}, {"_11_Q2_2.0_2.5", {2.0, 2.5}}};
        xrangelims = {{"_x_0_0.25", {0.00001, 0.25}}};
        zrangelims = {{"_z_0.3_0.35", {0.3, 0.35}}, {"_z_0.35_0.40", {0.350, 0.40}},{"_z_0.40_0.45",{0.40,0.45}},{"_z_0.45_0.5", {0.45, 0.5}}, {"_z_0.5_0.55", {0.50, 0.55}},{"_z_0.55_0.6",{0.55,0.6}}};

        // Vs ( = sqrt(2*M*E) for fixed target experiments, here E = 11 GeV)
        Vs = 4.7;
        //const double Vs = 28.635642;
      }
      else if(tab=="numbers_x2_pip_11_He3.txt" or tab=="numbers_x2_pim_11_He3.txt"){
        //std::cout<<" check condition 2"<<std::endl;
        // Initialize naming map for the Q-integration ranges (the first element is for the name of the output data file)
        Q2rangelims = {{"_11_Q2_1.50_2.0", {1.50, 2.0}}, {"_11_Q2_2.0_2.5", {2.0, 2.5}},{"_11_Q2_2.5_3.0",{2.5,3.0}}};
        xrangelims = {{"_x_0.25_0.5", {0.25, 0.5}}};
        zrangelims = {{"_z_0.30_0.40", {0.30, 0.40}}, {"_z_0.40_0.50", {0.40, 0.50}},{"_z_0.50_0.6",{0.50,0.60}}};

        // Vs ( = sqrt(2*M*E) for fixed target experiments, here E = 11 GeV)
        Vs = 4.7;
        //const double Vs = 28.635642;
      }
      else if(tab=="numbers_hQ2_x2_pip_11_He3.txt" or tab=="numbers_hQ2_x2_pim_11_He3.txt"){
        //std::cout<<" check condition 3"<<std::endl;
        // Initialize naming map for the Q-integration ranges (the first element is for the name of the output data file)
        Q2rangelims = {{"_11_Q2_3.0_8.0", {3.0, 8.0}}};
        xrangelims = {{"_x_0.25_0.70", {0.25, 0.70}}};
        zrangelims = {{"_z_0.30_0.40", {0.30, 0.40}}, {"_z_0.40_0.60", {0.40, 0.60}}};

        // Vs ( = sqrt(2*M*E) for fixed target experiments, here E = 11 GeV)
        Vs = 4.7;
        //const double Vs = 28.635642;
      }
      else if(tab=="numbers_x1_pip_11_Q22.5-3_He3.txt" or tab == "numbers_x1_pim_11_Q22.5-3_He3.txt"){
        // Initialize naming map for the Q-integration ranges (the first element is for the name of the output data file)
        Q2rangelims = {{"_11_Q2_2.5_3.0", {2.5, 3.0}}};
        xrangelims = {{"_x_0_0.25", {0.0,0.25}}};
        zrangelims = {{"_z_0.30_0.40", {0.30, 0.40}}, {"_z_0.40_0.50", {0.40, 0.50}},{"_z_0.50_0.6",{0.50,0.60}}};

        // Vs ( = sqrt(2*M*E) for fixed target experiments, here E = 11 GeV)
        Vs = 4.7;
      }
      //std::cout<<" file name "<<tab<<std::endl;
      // Create directory
      std::string opath = ProcessedDataPath + "/" + ofolder;
      mkdir(opath.c_str(), ACCESSPERMS);
      
      // Reading table with YAML
      const YAML::Node exp = YAML::LoadFile(RawDataFolder + "/" + tab);

      // Create an input filestream
      std::ifstream rawfile(RawDataFolder + tab);

      // Make sure the file is open
      if(!rawfile.is_open()) throw std::runtime_error("Could not open file");

      std::string line;
      // Vector to hold the rows of data
      std::vector<DataRow> rows;
      if(rawfile.good())
      {
        while (std::getline(rawfile,line)){
          // If the line starts whith i, skip it.
          // (They are comments in the header of the rawdata file).
          if (line[0] == 'i')
            continue;
          rows.push_back(parseLine(line));
        }

      }
      

      // Conditions to separate the values in the columns into different files.
      // Loop on Q bin boundaries.
      for (auto const& Q2 : Q2rangelims)
      {
        //std::cout<<Q2.first<<"check"<<std::endl;
        for (auto const& x : xrangelims)
        {
          for (auto const& z : zrangelims)
          {
            //std::cout<<" check numbers x "<<x.second.first<<"  "<< x.second.second <<" z "<<z.second.first <<"  "<< z.second.second <<" Q2 "<< Q2.second.first <<"   "<< Q2.second.second<<std::endl;
            int charge =0;
            std::string name;
            if (rows[0].had=="pi+") {charge=1;name="pip";}
            else if (rows[0].had=="pi-") {charge=-1;name="pim";}
            
            // File name picks up the m (i.e. Q) boundaries
            std::string ofileQ = "SoLID_Unpol" + Q2.first + x.first + z.first+"_"+name;

            // Initialize indexes vector for future selection
            // Vector to hold the indices of the selected rows
            std::vector<int> indexes;
            // Filter the rows based on the ranges
            for (const auto& row : rows) {
              //std::cout<<" check numbers"<<x.second.first<<" "<< row.x <<" "<< x.second.second <<" "<<z.second.first <<" "<<  row.z <<" "<< z.second.second <<" "<< Q2.second.first <<" "<< row.Q2 <<" "<< row.Q2 <<" "<< Q2.second.second<<std::endl;
              if (x.second.first < row.x && row.x < x.second.second && z.second.first < row.z && row.z < z.second.second && Q2.second.first < row.Q2 && row.Q2 < Q2.second.second) {
                indexes.push_back(row.index);
                //std::cout<<"check index"<<row.index<<std::endl;
              }
            }


            // Initialize result maps
            std::map<int, double> fdcross, fdstat, fdpT,fx,fy,fz,fQ2;

            // Prepare (outer) map for each output data file
            std::map<std::string, std::map<int, double>> filedata;

            // Finally select cross sections and their uncertainties.
            for (int i : indexes)
            {
              fdcross[i] = rows[i].value;
              fdstat[i]  = rows[i].stat;
              fdpT[i]    = rows[i].pT;
              fx[i]      = rows[i].x;
              fy[i]      = rows[i].y;
              fz[i]      = rows[i].z;
              fQ2[i]     = rows[i].Q2;
              //std::cout<<"check i "<<i<<" cross "<<fdcross[i]<<"pT"<<fdpT[i]<<std::endl;
            }

            filedata["cross"] = fdcross;
            filedata["stat"]  = fdstat;
            filedata["pT"]    = fdpT;
            filedata["x"]    = fx;
            filedata["y"]    = fy;
            filedata["z"]    = fz;
            filedata["Q2"]    = fQ2;

            // Plot labels
            std::map<std::string, std::string> labels
            {
              //{"xlabel", "#it{q}_{T} [GeV]"},
              //{"ylabel", "#frac{d^2#it{#sigma}}{d#it{Q}#d#it{q}_{T}} [pb GeV^{-2}]"},
              {"title", "SoLID, 3He target " + std::to_string(Q2.second.first) + " < Q2 < " + std::to_string(Q2.second.second)+","+ std::to_string(x.second.first) + " < x < " + std::to_string(x.second.second)+","+std::to_string(z.second.first) + " < z < " + std::to_string(z.second.second)},
                {"xlabelpy", "$P_{hT} \\rm{[GeV]}$"},
                {"ylabelpy", "$F^{\\pi^"+std::string(1,rows[0].had[2])+"}_{\\rm UU,T }\\left(x, z, q_T, Q^2 \\right)$"},
                {"titlepy", "SoLID unpol, " + std::to_string(Q2.second.first) + " < Q2 < " + std::to_string(Q2.second.second)}
            };

            /*
               NOTE on the conversion factor for the cross section.
               The raw data have a cross section expressed in cm**2/GeV**2/nucleon, but we would like to have it in NB./NUCLEON/GEV**2,
               then we have to convert: 1barn = 10**{-28}m**2 = 10**{-24}cm**2. Therefore, 1 cm**2 = 10**{24}barn= 10**{33}nb.

               NOTE on the calculation of y_min and y_max: y = arcsinh(sqrt{s}*xF/(2Q)).
               The value of x_min = 0, then y_min=0 for all bin in Q
               x_max = 1, then y_max = arcsinh(sqrt(s)/(2 Q_min)) for a specific Qmin < Q < Qmax bin
               */

            // Allocate emitter
            YAML::Emitter emit;

            // Write kinematics on the YAML emitter
            emit.SetFloatPrecision(8);
            emit.SetDoublePrecision(8);
            emit << YAML::BeginMap;
            emit << YAML::Key << "dependent_variables";
            emit << YAML::BeginSeq;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Value << YAML::Flow << labels;
            emit << YAML::Key << "qualifiers" << YAML::Value;
            emit << YAML::BeginSeq;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "process" << YAML::Key << "value" << YAML::Value << "SIDIS" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "observable" << YAML::Key << "value" << YAML::Value << "FUUT" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "target_isoscalarity" << YAML::Key << "value" << YAML::Value << "2/3" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "hadron" << YAML::Key << "value" << YAML::Value << "PI" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "charge" << YAML::Key << "value" << YAML::Value << charge << YAML::EndMap;
            //emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "beam" << YAML::Key << "value" << YAML::Value << rows[0].Eb << YAML::EndMap;
            //emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "hadron" << YAML::Key << "value" << YAML::Value << rows[0].had << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "prefactor" << YAML::Key << "value" << YAML::Value << 1 << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "Vs" << YAML::Key << "value" << YAML::Value << Vs << YAML::EndMap;
            //emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "Q2" << YAML::Key << "low" << YAML::Value << Q2.second.first << YAML::Key << "high" << YAML::Value << Q2.second.second  << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "Q" << YAML::Key << "low" << YAML::Value << sqrt(Q2.second.first) << YAML::Key << "high" << YAML::Value << sqrt(Q2.second.second)  << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            //emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "y" << YAML::Key
            //  << "low" << YAML::Value << asinh(0) << YAML::Key << "high" << YAML::Value << asinh(Vs / (2 * Q2.second.first)) << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "x" << YAML::Key
              << "low" << YAML::Value << x.second.first << YAML::Key << "high" << YAML::Value << x.second.second  << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "z" << YAML::Key
              << "low" << YAML::Value << z.second.first << YAML::Key << "high" << YAML::Value << z.second.second  << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "PS_reduction" << YAML::Key
               << "W" << YAML::Value << "2.3" << YAML::Key << "ymin" << YAML::Value << "0.01" << YAML::Key << "ymax" << YAML::Value << "0.95" << YAML::EndMap; 
               //<< "pTmin" << YAML::Value << "###" << YAML::Key << "etamin" << YAML::Value << "###" << YAML::Key << "etamax" << YAML::Value << "###" << YAML::EndMap; 
            emit << YAML::EndSeq;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& m : filedata["cross"])
            {
              //std::cout<<"check m"<<m.first<<std::endl;
              emit << YAML::BeginMap << YAML::Key << "errors" << YAML::Value << YAML::BeginSeq;

              // Oss: the conversion from cm**2 to barn is: 1b = 10**{-24}cm**2,
              //      but we need the cross section in pb: 1pb = 10**{-36}cm**2,
              //      <=> 1 cm**2 = 10**{36}pb
              // Remember conversion factor for the cross section
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "label" << YAML::Value << "unc" << YAML::Key << "value" << YAML::Value << (filedata["stat"][m.first])  << YAML::EndMap;
              //std::cout<<"check" <<filedata["stat"][m.first]<<std::endl;

              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "label" << YAML::Value << "add" << YAML::Key << "value" << YAML::Value << 0.1 << YAML::EndMap;
              emit << YAML::EndSeq;
              emit << YAML::Key << "value" << YAML::Value << (m.second) ;
              // emit << YAML::Key << "id"    << YAML::Value << m.first;
              emit << YAML::EndMap;
            }
            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::Key << "independent_variables";
            emit << YAML::BeginSeq;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Value << YAML::Flow;
            emit << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "PT" << YAML::Key << "units" << YAML::Value << "GEV" << YAML::EndMap;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            // for (auto const& p : filedata["pT"])
            //   {
            //     emit << YAML::Flow << YAML::BeginMap << YAML::Key << "value" << YAML::Value << p.second << YAML::Key << "high" << YAML::Value << (p.second + 0.125) << YAML::Key << "low" << YAML::Value << (p.second - 0.125) << YAML::Key
            //     /* Since Nangaparbat cross section is differential in qT, while this experimental set is differential in qT**2 and Q, we need to correct the value of the predictions by a factor DqT/(DQ * D(qT**2))*/
            //     << "factor" << YAML::Value << (p.second - p.first) / ((Q2.second.second - Q2.second.first) * (pow(p.second, 2) - pow(p.first, 2)))<< YAML::EndMap;
            //     // emit << YAML::Flow << YAML::BeginMap << YAML::Key << "value" << YAML::Value << p.second << YAML::EndMap;
            //   }
            for (auto const& p : filedata["pT"])
            {
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "value" << YAML::Value << p.second << YAML::Key
                /* Since Nangaparbat cross section is differential in qT, while this experimental set is differential in qT**2 and Q, we need to correct the value of the predictions by a factor DqT/(DQ * D(qT**2))*/
                << "factor" << YAML::Value << 1.0 / ((2*p.second)*(sqrt((Q2.second.second - Q2.second.first)/2)/((x.second.second-x.second.first)/2)/rows[0].Eb/0.93827208816))<< YAML::EndMap;
              //std::cout<<"check p.second"<<p.second<<std::endl;
              // emit << YAML::Flow << YAML::BeginMap << YAML::Key << "value" << YAML::Value << p.second << YAML::EndMap;
            }

            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "x" << YAML::EndMap;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& ix : filedata["x"])
              //std::cout<<" check x "<<ix.second<<std::endl;
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "high" << YAML::Value << x.second.second << YAML::Key << "low" << YAML::Value << x.second.first << YAML::Key << "value" << YAML::Value << ix.second << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "z" << YAML::EndMap;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& iz : filedata["z"])
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "high" << YAML::Value << z.second.second << YAML::Key << "low" << YAML::Value << z.second.first << YAML::Key << "value" << YAML::Value << iz.second << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "y" << YAML::EndMap;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& iy : filedata["y"])
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "high" << YAML::Value << 0.95 << YAML::Key << "low" << YAML::Value << 0.01 << YAML::Key << "value" << YAML::Value << iy.second << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "Q2" << YAML::EndMap;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& iQ2 : filedata["Q2"])
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "high" << YAML::Value << Q2.second.second << YAML::Key << "low" << YAML::Value << Q2.second.first << YAML::Key << "value" << YAML::Value << iQ2.second << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::EndMap;

            // Dump table to file
            std::ofstream fout(opath + "/" + ofileQ + ".yaml");
            fout << emit.c_str() << std::endl;
            fout.close();

            filenames.push_back(ofileQ);
          }

        }//loop over Q2
      }//loop overr x
    }//loop over z
       
    // Map to space properly the dataset names in datasets.yaml
    std::map<int, std::string> spaces = {{51," "}, {50,"  "}, {49,"   "}, {48, "    "}, {47, "     "}, {46, "      "},{45, "       "},{44,"        "}};

    // Produce outputnames to put in datasets.yaml
    std::string outputnames;
    for (std::string name : filenames)
      outputnames += "  - {name: " + name + "," + spaces[name.size()] + "file: " + name + ".yaml}\n";

    return outputnames;


    //return
    //  "check";
  }
}
namespace NangaParbat
{

  //_________________________________________________________________________________
  std::string PreprocessSoLIDKaon(std::string const& RawDataPath, std::string const& ProcessedDataPath)
  {
    std::cout << "Processing SoLID Kaon projection ..." << std::endl;

    // Path to the raw-data folder
    const std::string RawDataFolder = RawDataPath + "/SoLID_Unpol/";


    // Vector of tables to process
    const std::vector<std::string> tables = {"numbers_x1_pip_11_He3_kaon.txt","numbers_x1_pim_11_He3_kaon.txt","numbers_x2_pip_11_He3_kaon.txt","numbers_x2_pim_11_He3_kaon.txt","numbers_hQ2_x2_pip_11_He3_kaon.txt","numbers_hQ2_x2_pim_11_He3_kaon.txt","numbers_x1_pip_11_Q22.5-3_He3_kaon.txt","numbers_x1_pim_11_Q22.5-3_He3_kaon.txt"};

    // Output folder
    const std::string ofolder = "oSoLID_Unpol";

    // Vector for output names
    std::vector<std::string> filenames;

    // Loop over tables
    for (auto const& tab : tables)
    {
      //naming map
      std::map<std::string, std::pair<double, double>> Q2rangelims;
      std::map<std::string, std::pair<double, double>> xrangelims;
      std::map<std::string, std::pair<double, double>> zrangelims;
      
      //std::cout<<"check names "<<tab<<std::endl;
      // Vs ( = sqrt(2*M*E) for fixed target experiments)
      double Vs;
      if(tab=="numbers_x1_pip_11_He3_kaon.txt" or tab=="numbers_x1_pim_11_He3_kaon.txt"){
        //std::cout<<" check condition 1"<<std::endl;
        // Initialize naming map for the Q-integration ranges (the first element is for the name of the output data file)
        Q2rangelims = {{"_11_Q2_1_1.50", {1.0, 1.50}}, {"_11_Q2_1.50_2.0", {1.50, 2.0}}, {"_11_Q2_2.0_2.5", {2.0, 2.5}}};
        xrangelims = {{"_x_0_0.25", {0.00001, 0.25}}};
        zrangelims = {{"_z_0.3_0.35", {0.3, 0.35}}, {"_z_0.35_0.40", {0.350, 0.40}},{"_z_0.40_0.45",{0.40,0.45}},{"_z_0.45_0.5", {0.45, 0.5}}, {"_z_0.5_0.55", {0.50, 0.55}},{"_z_0.55_0.6",{0.55,0.6}}};

        // Vs ( = sqrt(2*M*E) for fixed target experiments, here E = 11 GeV)
        Vs = 4.7;
        //const double Vs = 28.635642;
      }
      else if(tab=="numbers_x2_pip_11_He3_kaon.txt" or tab=="numbers_x2_pim_11_He3_kaon.txt"){
        //std::cout<<" check condition 2"<<std::endl;
        // Initialize naming map for the Q-integration ranges (the first element is for the name of the output data file)
        Q2rangelims = {{"_11_Q2_1.50_2.0", {1.50, 2.0}}, {"_11_Q2_2.0_2.5", {2.0, 2.5}},{"_11_Q2_2.5_3.0",{2.5,3.0}}};
        xrangelims = {{"_x_0.25_0.5", {0.25, 0.5}}};
        zrangelims = {{"_z_0.30_0.40", {0.30, 0.40}}, {"_z_0.40_0.50", {0.40, 0.50}},{"_z_0.50_0.6",{0.50,0.60}}};

        // Vs ( = sqrt(2*M*E) for fixed target experiments, here E = 11 GeV)
        Vs = 4.7;
        //const double Vs = 28.635642;
      }
      else if(tab=="numbers_hQ2_x2_pip_11_He3_kaon.txt" or tab=="numbers_hQ2_x2_pim_11_He3_kaon.txt"){
        //std::cout<<" check condition 3"<<std::endl;
        // Initialize naming map for the Q-integration ranges (the first element is for the name of the output data file)
        Q2rangelims = {{"_11_Q2_3.0_8.0", {3.0, 8.0}}};
        xrangelims = {{"_x_0.25_0.70", {0.25, 0.70}}};
        zrangelims = {{"_z_0.30_0.40", {0.30, 0.40}}, {"_z_0.40_0.60", {0.40, 0.60}}};

        // Vs ( = sqrt(2*M*E) for fixed target experiments, here E = 11 GeV)
        Vs = 4.7;
        //const double Vs = 28.635642;
      }
      else if(tab=="numbers_x1_pip_11_Q22.5-3_He3_kaon.txt" or tab=="numbers_x1_pim_11_Q22.5-3_He3_kaon.txt"){
        // Initialize naming map for the Q-integration ranges (the first element is for the name of the output data file)
        Q2rangelims = {{"_11_Q2_2.5_3.0", {2.5, 3.0}}};
        xrangelims = {{"_x_0_0.25", {0.0,0.25}}};
        zrangelims = {{"_z_0.30_0.40", {0.30, 0.40}}, {"_z_0.40_0.50", {0.40, 0.50}},{"_z_0.50_0.6",{0.50,0.60}}};

        // Vs ( = sqrt(2*M*E) for fixed target experiments, here E = 11 GeV)
        Vs = 4.7;
      }
      // Create directory
      std::string opath = ProcessedDataPath + "/" + ofolder;
      mkdir(opath.c_str(), ACCESSPERMS);
      
      // Reading table with YAML
      const YAML::Node exp = YAML::LoadFile(RawDataFolder + "/" + tab);

      // Create an input filestream
      std::ifstream rawfile(RawDataFolder + tab);

      // Make sure the file is open
      if(!rawfile.is_open()) throw std::runtime_error("Could not open file");

      std::string line;
      // Vector to hold the rows of data
      std::vector<DataRow> rows;
      if(rawfile.good())
      {
        while (std::getline(rawfile,line)){
          // If the line starts whith i, skip it.
          // (They are comments in the header of the rawdata file).
          if (line[0] == 'i')
            continue;
          rows.push_back(parseLine(line));
        }

      }
      

      // Conditions to separate the values in the columns into different files.
      // Loop on Q bin boundaries.
      for (auto const& Q2 : Q2rangelims)
      {
        //std::cout<<Q2.first<<"check"<<std::endl;
        for (auto const& x : xrangelims)
        {
          for (auto const& z : zrangelims)
          {
            //std::cout<<" check numbers"<<x.second.first<<"  "<< x.second.second <<" "<<z.second.first <<"  "<< z.second.second <<" "<< Q2.second.first <<"   "<< Q2.second.second<<std::endl;
            int charge =0;
            std::string name;
            if (rows[0].had=="K+") {charge=1;name="Kp";}
            else if (rows[0].had=="K-") {charge=-1;name="Km";}
            
            // File name picks up the m (i.e. Q) boundaries
            std::string ofileQ = "SoLID_Unpol" + Q2.first + x.first + z.first+"_"+name;

            // Initialize indexes vector for future selection
            // Vector to hold the indices of the selected rows
            std::vector<int> indexes;
            // Filter the rows based on the ranges
            for (const auto& row : rows) {
              //std::cout<<" check numbers"<<x.second.first<<" "<< row.x <<" "<< x.second.second <<" "<<z.second.first <<" "<<  row.z <<" "<< z.second.second <<" "<< Q2.second.first <<" "<< row.Q2 <<" "<< row.Q2 <<" "<< Q2.second.second<<std::endl;
              if (x.second.first < row.x && row.x < x.second.second && z.second.first < row.z && row.z < z.second.second && Q2.second.first < row.Q2 && row.Q2 < Q2.second.second) {
                indexes.push_back(row.index);
                //std::cout<<"check index"<<row.index<<std::endl;
              }
            }


            // Initialize result maps
            std::map<int, double> fdcross, fdstat, fdpT,fx,fy,fz,fQ2;

            // Prepare (outer) map for each output data file
            std::map<std::string, std::map<int, double>> filedata;

            // Finally select cross sections and their uncertainties.
            for (int i : indexes)
            {
              fdcross[i] = rows[i].value;
              fdstat[i]  = rows[i].stat;
              fdpT[i]    = rows[i].pT;
              fx[i]      = rows[i].x;
              fy[i]      = rows[i].y;
              fz[i]      = rows[i].z;
              fQ2[i]     = rows[i].Q2;
              //std::cout<<"check i "<<i<<" cross "<<fdcross[i]<<"pT"<<fdpT[i]<<std::endl;
            }

            filedata["cross"] = fdcross;
            filedata["stat"]  = fdstat;
            filedata["pT"]    = fdpT;
            filedata["x"]    = fx;
            filedata["y"]    = fy;
            filedata["z"]    = fz;
            filedata["Q2"]    = fQ2;

            // Plot labels
            std::map<std::string, std::string> labels
            {
              //{"xlabel", "#it{q}_{T} [GeV]"},
              //{"ylabel", "#frac{d^2#it{#sigma}}{d#it{Q}#d#it{q}_{T}} [pb GeV^{-2}]"},
              {"title", "SoLID, 3He target " + std::to_string(Q2.second.first) + " < Q2 < " + std::to_string(Q2.second.second)+","+ std::to_string(x.second.first) + " < x < " + std::to_string(x.second.second)+","+std::to_string(z.second.first) + " < z < " + std::to_string(z.second.second)},
                {"xlabelpy", "$P_{hT} \\rm{[GeV]}$"},
                {"ylabelpy", "$F^{K^"+std::string(1,rows[0].had[1])+"}_{\\rm UU,T }\\left(x, z, q_T, Q^2 \\right)$"},
                {"titlepy", "SoLID unpol, " + std::to_string(Q2.second.first) + " < Q2 < " + std::to_string(Q2.second.second)}
            };

            /*
               NOTE on the conversion factor for the cross section.
               The raw data have a cross section expressed in cm**2/GeV**2/nucleon, but we would like to have it in NB./NUCLEON/GEV**2,
               then we have to convert: 1barn = 10**{-28}m**2 = 10**{-24}cm**2. Therefore, 1 cm**2 = 10**{24}barn= 10**{33}nb.

               NOTE on the calculation of y_min and y_max: y = arcsinh(sqrt{s}*xF/(2Q)).
               The value of x_min = 0, then y_min=0 for all bin in Q
               x_max = 1, then y_max = arcsinh(sqrt(s)/(2 Q_min)) for a specific Qmin < Q < Qmax bin
               */

            // Allocate emitter
            YAML::Emitter emit;

            // Write kinematics on the YAML emitter
            emit.SetFloatPrecision(8);
            emit.SetDoublePrecision(8);
            emit << YAML::BeginMap;
            emit << YAML::Key << "dependent_variables";
            emit << YAML::BeginSeq;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Value << YAML::Flow << labels;
            emit << YAML::Key << "qualifiers" << YAML::Value;
            emit << YAML::BeginSeq;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "process" << YAML::Key << "value" << YAML::Value << "SIDIS" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "observable" << YAML::Key << "value" << YAML::Value << "FUUT" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "target_isoscalarity" << YAML::Key << "value" << YAML::Value << "2/3" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "hadron" << YAML::Key << "value" << YAML::Value << "KAON" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "charge" << YAML::Key << "value" << YAML::Value << charge << YAML::EndMap;
            //emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "beam" << YAML::Key << "value" << YAML::Value << rows[0].Eb << YAML::EndMap;
            //emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "hadron" << YAML::Key << "value" << YAML::Value << rows[0].had << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "prefactor" << YAML::Key << "value" << YAML::Value << 1 << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "Vs" << YAML::Key << "value" << YAML::Value << Vs << YAML::EndMap;
            //emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "Q2" << YAML::Key << "low" << YAML::Value << Q2.second.first << YAML::Key << "high" << YAML::Value << Q2.second.second  << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "Q" << YAML::Key << "low" << YAML::Value << sqrt(Q2.second.first) << YAML::Key << "high" << YAML::Value << sqrt(Q2.second.second)  << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            //emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "y" << YAML::Key
            //  << "low" << YAML::Value << asinh(0) << YAML::Key << "high" << YAML::Value << asinh(Vs / (2 * Q2.second.first)) << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "x" << YAML::Key
              << "low" << YAML::Value << x.second.first << YAML::Key << "high" << YAML::Value << x.second.second  << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "z" << YAML::Key
              << "low" << YAML::Value << z.second.first << YAML::Key << "high" << YAML::Value << z.second.second  << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "PS_reduction" << YAML::Key
               << "W" << YAML::Value << "2.3" << YAML::Key << "ymin" << YAML::Value << "0.01" << YAML::Key << "ymax" << YAML::Value << "0.95" << YAML::EndMap; 
               //<< "pTmin" << YAML::Value << "###" << YAML::Key << "etamin" << YAML::Value << "###" << YAML::Key << "etamax" << YAML::Value << "###" << YAML::EndMap; 
            emit << YAML::EndSeq;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& m : filedata["cross"])
            {
              //std::cout<<"check m"<<m.first<<std::endl;
              emit << YAML::BeginMap << YAML::Key << "errors" << YAML::Value << YAML::BeginSeq;

              // Oss: the conversion from cm**2 to barn is: 1b = 10**{-24}cm**2,
              //      but we need the cross section in pb: 1pb = 10**{-36}cm**2,
              //      <=> 1 cm**2 = 10**{36}pb
              // Remember conversion factor for the cross section
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "label" << YAML::Value << "unc" << YAML::Key << "value" << YAML::Value << (filedata["stat"][m.first])  << YAML::EndMap;
              //std::cout<<"check" <<filedata["stat"][m.first]<<std::endl;

              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "label" << YAML::Value << "add" << YAML::Key << "value" << YAML::Value << 0.1 << YAML::EndMap;
              emit << YAML::EndSeq;
              emit << YAML::Key << "value" << YAML::Value << (m.second) ;
              // emit << YAML::Key << "id"    << YAML::Value << m.first;
              emit << YAML::EndMap;
            }
            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::Key << "independent_variables";
            emit << YAML::BeginSeq;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Value << YAML::Flow;
            emit << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "PT" << YAML::Key << "units" << YAML::Value << "GEV" << YAML::EndMap;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            // for (auto const& p : filedata["pT"])
            //   {
            //     emit << YAML::Flow << YAML::BeginMap << YAML::Key << "value" << YAML::Value << p.second << YAML::Key << "high" << YAML::Value << (p.second + 0.125) << YAML::Key << "low" << YAML::Value << (p.second - 0.125) << YAML::Key
            //     /* Since Nangaparbat cross section is differential in qT, while this experimental set is differential in qT**2 and Q, we need to correct the value of the predictions by a factor DqT/(DQ * D(qT**2))*/
            //     << "factor" << YAML::Value << (p.second - p.first) / ((Q2.second.second - Q2.second.first) * (pow(p.second, 2) - pow(p.first, 2)))<< YAML::EndMap;
            //     // emit << YAML::Flow << YAML::BeginMap << YAML::Key << "value" << YAML::Value << p.second << YAML::EndMap;
            //   }
            for (auto const& p : filedata["pT"])
            {
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "value" << YAML::Value << p.second << YAML::Key
                /* Since Nangaparbat cross section is differential in qT, while this experimental set is differential in qT**2 and Q, we need to correct the value of the predictions by a factor DqT/(DQ * D(qT**2))*/
                << "factor" << YAML::Value << 1.0 / ((2*p.second)*(sqrt((Q2.second.second - Q2.second.first)/2)/((x.second.second-x.second.first)/2)/rows[0].Eb/0.93827208816))<< YAML::EndMap;
              //std::cout<<"check p.second"<<p.second<<std::endl;
              // emit << YAML::Flow << YAML::BeginMap << YAML::Key << "value" << YAML::Value << p.second << YAML::EndMap;
            }

            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "x" << YAML::EndMap;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& ix : filedata["x"])
              //std::cout<<" check x "<<ix.second<<std::endl;
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "high" << YAML::Value << x.second.second << YAML::Key << "low" << YAML::Value << x.second.first << YAML::Key << "value" << YAML::Value << ix.second << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "z" << YAML::EndMap;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& iz : filedata["z"])
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "high" << YAML::Value << z.second.second << YAML::Key << "low" << YAML::Value << z.second.first << YAML::Key << "value" << YAML::Value << iz.second << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "y" << YAML::EndMap;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& iy : filedata["y"])
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "high" << YAML::Value << 0.95 << YAML::Key << "low" << YAML::Value << 0.01 << YAML::Key << "value" << YAML::Value << iy.second << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "Q2" << YAML::EndMap;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& iQ2 : filedata["Q2"])
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "high" << YAML::Value << Q2.second.second << YAML::Key << "low" << YAML::Value << Q2.second.first << YAML::Key << "value" << YAML::Value << iQ2.second << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::EndMap;

            // Dump table to file
            std::ofstream fout(opath + "/" + ofileQ + ".yaml");
            fout << emit.c_str() << std::endl;
            fout.close();

            filenames.push_back(ofileQ);
          }

        }//loop over Q2
      }//loop overr x
    }//loop over z
       
    // Map to space properly the dataset names in datasets.yaml
    std::map<int, std::string> spaces = {{51," "}, {50,"  "}, {49,"   "}, {48, "    "}, {47, "     "}, {46, "      "},{45, "       "},{44,"        "}};

    // Produce outputnames to put in datasets.yaml
    std::string outputnames;
    for (std::string name : filenames)
      outputnames += "  - {name: " + name + "," + spaces[name.size()] + "file: " + name + ".yaml}\n";

    return outputnames;


    //return
    //  "check";
  }
}
namespace NangaParbat
{

  //_________________________________________________________________________________
  std::string PreprocessSoLID8p8(std::string const& RawDataPath, std::string const& ProcessedDataPath)
  {
    std::cout << "Processing SoLID Pion projection for 8.8 GeV beam ..." << std::endl;

    // Path to the raw-data folder
    const std::string RawDataFolder = RawDataPath + "/SoLID_Unpol/";


    // Vector of tables to process
    const std::vector<std::string> tables = {"numbers_x1_pip_8p8_He3.txt","numbers_x1_pim_8p8_He3.txt","numbers_x2_pip_8p8_He3.txt","numbers_x2_pim_8p8_He3.txt","numbers_x2_pip_8p8_hQ2_He3.txt","numbers_x2_pim_8p8_hQ2_He3.txt","numbers_x1_pip_8p8_Q22-2.5_He3.txt","numbers_x1_pim_8p8_Q21-1.5_He3.txt"};

    // Output folder
    const std::string ofolder = "oSoLID_Unpol";

    // Vector for output names
    std::vector<std::string> filenames;

    // Loop over tables
    for (auto const& tab : tables)
    {
      //naming map
      std::map<std::string, std::pair<double, double>> Q2rangelims;
      std::map<std::string, std::pair<double, double>> xrangelims;
      std::map<std::string, std::pair<double, double>> zrangelims;
      
      //std::cout<<"check names "<<tab<<std::endl;
      // Vs ( = sqrt(2*M*E) for fixed target experiments)
      double Vs;
      if(tab=="numbers_x1_pip_8p8_He3.txt"){
        //std::cout<<" check condition 1"<<std::endl;
        // Initialize naming map for the Q-integration ranges (the first element is for the name of the output data file)
        Q2rangelims = {{"_8p8_Q2_1.0_1.5", {1.0, 1.5}}, {"_8p8_Q2_1.5_2", {1.5, 2}}};
        xrangelims = {{"_x_0_0.25", {0.00001, 0.25}}};
        zrangelims = {{"_z_0.3_0.35", {0.3, 0.35}}, {"_z_0.35_0.40", {0.350, 0.40}},{"_z_0.40_0.45",{0.40,0.45}},{"_z_0.45_0.5", {0.45, 0.5}}, {"_z_0.5_0.55", {0.50, 0.55}},{"_z_0.55_0.6",{0.55,0.6}}};

        // Vs ( = sqrt(2*M*E) for fixed target experiments, here E = 8p8 GeV)
        Vs = 4.0631;
      }
      else if(tab=="numbers_x1_pip_8p8_Q22-2.5_He3.txt"){
        Q2rangelims = {{"_8p8_Q2_2_2.5", {2, 2.5}}};
        xrangelims = {{"_x_0_0.25", {0.00001, 0.25}}};
        zrangelims = {{"_z_0.3_0.4", {0.30, 0.40}}, {"_z_0.4_0.60", {0.40, 0.60}}};

        // Vs ( = sqrt(2*M*E) for fixed target experiments, here E = 8p8 GeV)
        Vs = 4.0631;
      }
      else if(tab=="numbers_x2_pip_8p8_He3.txt" or tab=="numbers_x2_pim_8p8_He3.txt"){
        Q2rangelims = {{"_8p8_Q2_1.5_2", {1.5, 2}}, {"_8p8_Q2_2_2.5", {2, 2.5}},{"_8p8_Q2_2.5-3",{2.5,3}}};
        xrangelims = {{"_x_0.25_0.5", {0.25, 0.5}}};
        zrangelims = {{"_z_0.3_0.35", {0.3, 0.35}}, {"_z_0.35_0.40", {0.350, 0.40}},{"_z_0.40_0.45",{0.40,0.45}},{"_z_0.45_0.5", {0.45, 0.5}}, {"_z_0.5_0.55", {0.50, 0.55}},{"_z_0.55_0.6",{0.55,0.6}}};

        // Vs ( = sqrt(2*M*E) for fixed target experiments, here E = 8p8 GeV)
        Vs = 4.0631;
      }
      else if(tab=="numbers_x1_pim_8p8_Q21-1.5_He3.txt"){
        Q2rangelims = {{"_8p8_Q2_1.0_1.5", {1.0, 1.5}}};
        xrangelims = {{"_x_0_0.25", {0.00001, 0.25}}};
        zrangelims = {{"_z_0.3_0.35", {0.3, 0.35}}, {"_z_0.35_0.40", {0.350, 0.40}},{"_z_0.40_0.45",{0.40,0.45}},{"_z_0.45_0.5", {0.45, 0.5}}, {"_z_0.5_0.55", {0.50, 0.55}},{"_z_0.55_0.6",{0.55,0.6}}};
      }
      else if(tab=="numbers_x1_pim_8p8_He3.txt"){
        Q2rangelims = {{"_8p8_Q2_1.5_2", {1.5, 2}},{"_8p8_Q2_2_2.5", {2, 2.5}}};
        xrangelims = {{"_x_0_0.25", {0.00001, 0.25}}};
        zrangelims = {{"_z_0.3_0.4", {0.3, 0.4}}, {"_z_0.4_0.60", {0.40, 0.60}}};
      }
      else if(tab=="numbers_x2_pim_8p8_hQ2_He3.txt" or tab=="numbers_x2_pip_8p8_hQ2_He3.txt"){
        Q2rangelims = {{"_8p8_Q2_3_8", {3, 8}}};
        xrangelims = {{"_x_0.25_0.7", {0.25,07}}};
        zrangelims = {{"_z_0.3_0.4", {0.3, 0.4}}, {"_z_0.4_0.60", {0.40, 0.60}}};
      }
      //std::cout<<" file name "<<tab<<std::endl;
      // Create directory
      std::string opath = ProcessedDataPath + "/" + ofolder;
      mkdir(opath.c_str(), ACCESSPERMS);
      
      // Reading table with YAML
      const YAML::Node exp = YAML::LoadFile(RawDataFolder + "/" + tab);

      // Create an input filestream
      std::ifstream rawfile(RawDataFolder + tab);

      // Make sure the file is open
      if(!rawfile.is_open()) throw std::runtime_error("Could not open file");

      std::string line;
      // Vector to hold the rows of data
      std::vector<DataRow> rows;
      if(rawfile.good())
      {
        while (std::getline(rawfile,line)){
          // If the line starts whith i, skip it.
          // (They are comments in the header of the rawdata file).
          if (line[0] == 'i')
            continue;
          rows.push_back(parseLine(line));
        }

      }
      

      // Conditions to separate the values in the columns into different files.
      // Loop on Q bin boundaries.
      for (auto const& Q2 : Q2rangelims)
      {
        //std::cout<<Q2.first<<"check"<<std::endl;
        for (auto const& x : xrangelims)
        {
          for (auto const& z : zrangelims)
          {
            //std::cout<<" check numbers x "<<x.second.first<<"  "<< x.second.second <<" z "<<z.second.first <<"  "<< z.second.second <<" Q2 "<< Q2.second.first <<"   "<< Q2.second.second<<std::endl;
            int charge =0;
            std::string name;
            if (rows[0].had=="pi+") {charge=1;name="pip";}
            else if (rows[0].had=="pi-") {charge=-1;name="pim";}
            
            // File name picks up the m (i.e. Q) boundaries
            std::string ofileQ = "SoLID_Unpol" + Q2.first + x.first + z.first+"_"+name;

            // Initialize indexes vector for future selection
            // Vector to hold the indices of the selected rows
            std::vector<int> indexes;
            // Filter the rows based on the ranges
            for (const auto& row : rows) {
              //std::cout<<" check numbers"<<x.second.first<<" "<< row.x <<" "<< x.second.second <<" "<<z.second.first <<" "<<  row.z <<" "<< z.second.second <<" "<< Q2.second.first <<" "<< row.Q2 <<" "<< row.Q2 <<" "<< Q2.second.second<<std::endl;
              if (x.second.first < row.x && row.x < x.second.second && z.second.first < row.z && row.z < z.second.second && Q2.second.first < row.Q2 && row.Q2 < Q2.second.second) {
                indexes.push_back(row.index);
                //std::cout<<"check index"<<row.index<<std::endl;
              }
            }


            // Initialize result maps
            std::map<int, double> fdcross, fdstat, fdpT,fx,fy,fz,fQ2;

            // Prepare (outer) map for each output data file
            std::map<std::string, std::map<int, double>> filedata;

            // Finally select cross sections and their uncertainties.
            for (int i : indexes)
            {
              fdcross[i] = rows[i].value;
              fdstat[i]  = rows[i].stat;
              fdpT[i]    = rows[i].pT;
              fx[i]      = rows[i].x;
              fy[i]      = rows[i].y;
              fz[i]      = rows[i].z;
              fQ2[i]     = rows[i].Q2;
              //std::cout<<"check i "<<i<<" cross "<<fdcross[i]<<"pT"<<fdpT[i]<<std::endl;
            }

            filedata["cross"] = fdcross;
            filedata["stat"]  = fdstat;
            filedata["pT"]    = fdpT;
            filedata["x"]    = fx;
            filedata["y"]    = fy;
            filedata["z"]    = fz;
            filedata["Q2"]    = fQ2;

            // Plot labels
            std::map<std::string, std::string> labels
            {
              //{"xlabel", "#it{q}_{T} [GeV]"},
              //{"ylabel", "#frac{d^2#it{#sigma}}{d#it{Q}#d#it{q}_{T}} [pb GeV^{-2}]"},
              {"title", "SoLID, 3He target " + std::to_string(Q2.second.first) + " < Q2 < " + std::to_string(Q2.second.second)+","+ std::to_string(x.second.first) + " < x < " + std::to_string(x.second.second)+","+std::to_string(z.second.first) + " < z < " + std::to_string(z.second.second)},
                {"xlabelpy", "$P_{hT} \\rm{[GeV]}$"},
                {"ylabelpy", "$F^{\\pi^"+std::string(1,rows[0].had[2])+"}_{\\rm UU,T }\\left(x, z, q_T, Q^2 \\right)$"},
                {"titlepy", "SoLID unpol, " + std::to_string(Q2.second.first) + " < Q2 < " + std::to_string(Q2.second.second)}
            };

            /*
               NOTE on the conversion factor for the cross section.
               The raw data have a cross section expressed in cm**2/GeV**2/nucleon, but we would like to have it in NB./NUCLEON/GEV**2,
               then we have to convert: 1barn = 10**{-28}m**2 = 10**{-24}cm**2. Therefore, 1 cm**2 = 10**{24}barn= 10**{33}nb.

               NOTE on the calculation of y_min and y_max: y = arcsinh(sqrt{s}*xF/(2Q)).
               The value of x_min = 0, then y_min=0 for all bin in Q
               x_max = 1, then y_max = arcsinh(sqrt(s)/(2 Q_min)) for a specific Qmin < Q < Qmax bin
               */

            // Allocate emitter
            YAML::Emitter emit;

            // Write kinematics on the YAML emitter
            emit.SetFloatPrecision(8);
            emit.SetDoublePrecision(8);
            emit << YAML::BeginMap;
            emit << YAML::Key << "dependent_variables";
            emit << YAML::BeginSeq;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Value << YAML::Flow << labels;
            emit << YAML::Key << "qualifiers" << YAML::Value;
            emit << YAML::BeginSeq;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "process" << YAML::Key << "value" << YAML::Value << "SIDIS" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "observable" << YAML::Key << "value" << YAML::Value << "FUUT" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "target_isoscalarity" << YAML::Key << "value" << YAML::Value << "2/3" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "hadron" << YAML::Key << "value" << YAML::Value << "PI" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "charge" << YAML::Key << "value" << YAML::Value << charge << YAML::EndMap;
            //emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "beam" << YAML::Key << "value" << YAML::Value << rows[0].Eb << YAML::EndMap;
            //emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "hadron" << YAML::Key << "value" << YAML::Value << rows[0].had << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "prefactor" << YAML::Key << "value" << YAML::Value << 1 << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "Vs" << YAML::Key << "value" << YAML::Value << Vs << YAML::EndMap;
            //emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "Q2" << YAML::Key << "low" << YAML::Value << Q2.second.first << YAML::Key << "high" << YAML::Value << Q2.second.second  << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "Q" << YAML::Key << "low" << YAML::Value << sqrt(Q2.second.first) << YAML::Key << "high" << YAML::Value << sqrt(Q2.second.second)  << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            //emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "y" << YAML::Key
            //  << "low" << YAML::Value << asinh(0) << YAML::Key << "high" << YAML::Value << asinh(Vs / (2 * Q2.second.first)) << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "x" << YAML::Key
              << "low" << YAML::Value << x.second.first << YAML::Key << "high" << YAML::Value << x.second.second  << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "z" << YAML::Key
              << "low" << YAML::Value << z.second.first << YAML::Key << "high" << YAML::Value << z.second.second  << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "PS_reduction" << YAML::Key
               << "W" << YAML::Value << "2.3" << YAML::Key << "ymin" << YAML::Value << "0.01" << YAML::Key << "ymax" << YAML::Value << "0.95" << YAML::EndMap; 
               //<< "pTmin" << YAML::Value << "###" << YAML::Key << "etamin" << YAML::Value << "###" << YAML::Key << "etamax" << YAML::Value << "###" << YAML::EndMap; 
            emit << YAML::EndSeq;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& m : filedata["cross"])
            {
              //std::cout<<"check m"<<m.first<<std::endl;
              emit << YAML::BeginMap << YAML::Key << "errors" << YAML::Value << YAML::BeginSeq;

              // Oss: the conversion from cm**2 to barn is: 1b = 10**{-24}cm**2,
              //      but we need the cross section in pb: 1pb = 10**{-36}cm**2,
              //      <=> 1 cm**2 = 10**{36}pb
              // Remember conversion factor for the cross section
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "label" << YAML::Value << "unc" << YAML::Key << "value" << YAML::Value << (filedata["stat"][m.first])  << YAML::EndMap;
              //std::cout<<"check" <<filedata["stat"][m.first]<<std::endl;

              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "label" << YAML::Value << "add" << YAML::Key << "value" << YAML::Value << 0.1 << YAML::EndMap;
              emit << YAML::EndSeq;
              emit << YAML::Key << "value" << YAML::Value << (m.second) ;
              // emit << YAML::Key << "id"    << YAML::Value << m.first;
              emit << YAML::EndMap;
            }
            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::Key << "independent_variables";
            emit << YAML::BeginSeq;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Value << YAML::Flow;
            emit << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "PT" << YAML::Key << "units" << YAML::Value << "GEV" << YAML::EndMap;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            // for (auto const& p : filedata["pT"])
            //   {
            //     emit << YAML::Flow << YAML::BeginMap << YAML::Key << "value" << YAML::Value << p.second << YAML::Key << "high" << YAML::Value << (p.second + 0.125) << YAML::Key << "low" << YAML::Value << (p.second - 0.125) << YAML::Key
            //     /* Since Nangaparbat cross section is differential in qT, while this experimental set is differential in qT**2 and Q, we need to correct the value of the predictions by a factor DqT/(DQ * D(qT**2))*/
            //     << "factor" << YAML::Value << (p.second - p.first) / ((Q2.second.second - Q2.second.first) * (pow(p.second, 2) - pow(p.first, 2)))<< YAML::EndMap;
            //     // emit << YAML::Flow << YAML::BeginMap << YAML::Key << "value" << YAML::Value << p.second << YAML::EndMap;
            //   }
            for (auto const& p : filedata["pT"])
            {
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "value" << YAML::Value << p.second << YAML::Key
                /* Since Nangaparbat cross section is differential in qT, while this experimental set is differential in qT**2 and Q, we need to correct the value of the predictions by a factor DqT/(DQ * D(qT**2))*/
                << "factor" << YAML::Value << 1.0 / ((2*p.second)*(sqrt((Q2.second.second - Q2.second.first)/2)/((x.second.second-x.second.first)/2)/rows[0].Eb/0.93827208816))<< YAML::EndMap;
              //std::cout<<"check p.second"<<p.second<<std::endl;
              // emit << YAML::Flow << YAML::BeginMap << YAML::Key << "value" << YAML::Value << p.second << YAML::EndMap;
            }

            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "x" << YAML::EndMap;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& ix : filedata["x"])
              //std::cout<<" check x "<<ix.second<<std::endl;
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "high" << YAML::Value << x.second.second << YAML::Key << "low" << YAML::Value << x.second.first << YAML::Key << "value" << YAML::Value << ix.second << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "z" << YAML::EndMap;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& iz : filedata["z"])
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "high" << YAML::Value << z.second.second << YAML::Key << "low" << YAML::Value << z.second.first << YAML::Key << "value" << YAML::Value << iz.second << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "y" << YAML::EndMap;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& iy : filedata["y"])
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "high" << YAML::Value << 0.95 << YAML::Key << "low" << YAML::Value << 0.01 << YAML::Key << "value" << YAML::Value << iy.second << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "Q2" << YAML::EndMap;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& iQ2 : filedata["Q2"])
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "high" << YAML::Value << Q2.second.second << YAML::Key << "low" << YAML::Value << Q2.second.first << YAML::Key << "value" << YAML::Value << iQ2.second << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::EndMap;

            // Dump table to file
            std::ofstream fout(opath + "/" + ofileQ + ".yaml");
            fout << emit.c_str() << std::endl;
            fout.close();

            filenames.push_back(ofileQ);
          }

        }//loop over Q2
      }//loop overr x
    }//loop over z
       
    // Map to space properly the dataset names in datasets.yaml
    std::map<int, std::string> spaces = {{51," "}, {50,"  "}, {49,"   "}, {48, "    "}, {47, "     "}, {46, "      "},{45, "       "},{44,"        "}};

    // Produce outputnames to put in datasets.yaml
    std::string outputnames;
    for (std::string name : filenames)
      outputnames += "  - {name: " + name + "," + spaces[name.size()] + "file: " + name + ".yaml}\n";

    return outputnames;


    //return
    //  "check";
  }
}
namespace NangaParbat
{

  //_________________________________________________________________________________
  std::string PreprocessSoLID8p8Kaon(std::string const& RawDataPath, std::string const& ProcessedDataPath)
  {
    std::cout << "Processing SoLID Kaon projection for 8.8 GeV beam ..." << std::endl;

    // Path to the raw-data folder
    const std::string RawDataFolder = RawDataPath + "/SoLID_Unpol/";


    // Vector of tables to process
    const std::vector<std::string> tables = {"numbers_x1_pip_8p8_He3_kaon.txt","numbers_x1_pim_8p8_He3_kaon.txt","numbers_x2_pip_8p8_He3_kaon.txt","numbers_x2_pim_8p8_He3_kaon.txt","numbers_x2_pip_8p8_hQ2_He3_kaon.txt","numbers_x2_pim_8p8_hQ2_He3_kaon.txt","numbers_x1_pip_8p8_Q22-2.5_He3_kaon.txt","numbers_x1_pim_8p8_Q21-1.5_He3_kaon.txt"};

    // Output folder
    const std::string ofolder = "oSoLID_Unpol";

    // Vector for output names
    std::vector<std::string> filenames;

    // Loop over tables
    for (auto const& tab : tables)
    {
      //naming map
      std::map<std::string, std::pair<double, double>> Q2rangelims;
      std::map<std::string, std::pair<double, double>> xrangelims;
      std::map<std::string, std::pair<double, double>> zrangelims;
      
      //std::cout<<"check names "<<tab<<std::endl;
      // Vs ( = sqrt(2*M*E) for fixed target experiments)
      double Vs;
      if(tab=="numbers_x1_pip_8p8_He3_kaon.txt"){
        //std::cout<<" check condition 1"<<std::endl;
        // Initialize naming map for the Q-integration ranges (the first element is for the name of the output data file)
        Q2rangelims = {{"_8p8_Q2_1.0_1.5", {1.0, 1.5}}, {"_8p8_Q2_1.5_2", {1.5, 2}}};
        xrangelims = {{"_x_0_0.25", {0.00001, 0.25}}};
        zrangelims = {{"_z_0.3_0.35", {0.3, 0.35}}, {"_z_0.35_0.40", {0.350, 0.40}},{"_z_0.40_0.45",{0.40,0.45}},{"_z_0.45_0.5", {0.45, 0.5}}, {"_z_0.5_0.55", {0.50, 0.55}},{"_z_0.55_0.6",{0.55,0.6}}};

        // Vs ( = sqrt(2*M*E) for fixed target experiments, here E = 8p8 GeV)
        Vs = 4.0631;
      }
      else if(tab=="numbers_x1_pip_8p8_Q22-2.5_He3_kaon.txt"){
        Q2rangelims = {{"_8p8_Q2_2_2.5", {2, 2.5}}};
        xrangelims = {{"_x_0_0.25", {0.00001, 0.25}}};
        zrangelims = {{"_z_0.3_0.4", {0.30, 0.40}}, {"_z_0.4_0.60", {0.40, 0.60}}};

        // Vs ( = sqrt(2*M*E) for fixed target experiments, here E = 8p8 GeV)
        Vs = 4.0631;
      }
      else if(tab=="numbers_x2_pip_8p8_He3_kaon.txt" or tab=="numbers_x2_pim_8p8_He3_kaon.txt"){
        Q2rangelims = {{"_8p8_Q2_1.5_2", {1.5, 2}}, {"_8p8_Q2_2_2.5", {2, 2.5}},{"_8p8_Q2_2.5-3",{2.5,3}}};
        xrangelims = {{"_x_0.25_0.5", {0.25, 0.5}}};
        zrangelims = {{"_z_0.3_0.35", {0.3, 0.35}}, {"_z_0.35_0.40", {0.350, 0.40}},{"_z_0.40_0.45",{0.40,0.45}},{"_z_0.45_0.5", {0.45, 0.5}}, {"_z_0.5_0.55", {0.50, 0.55}},{"_z_0.55_0.6",{0.55,0.6}}};

        // Vs ( = sqrt(2*M*E) for fixed target experiments, here E = 8p8 GeV)
        Vs = 4.0631;
      }
      else if(tab=="numbers_x1_pim_8p8_Q21-1.5_He3_kaon.txt"){
        Q2rangelims = {{"_8p8_Q2_1.0_1.5", {1.0, 1.5}}};
        xrangelims = {{"_x_0_0.25", {0.00001, 0.25}}};
        zrangelims = {{"_z_0.3_0.35", {0.3, 0.35}}, {"_z_0.35_0.40", {0.350, 0.40}},{"_z_0.40_0.45",{0.40,0.45}},{"_z_0.45_0.5", {0.45, 0.5}}, {"_z_0.5_0.55", {0.50, 0.55}},{"_z_0.55_0.6",{0.55,0.6}}};
      }
      else if(tab=="numbers_x1_pim_8p8_He3_kaon.txt"){
        Q2rangelims = {{"_8p8_Q2_1.5_2", {1.5, 2}},{"_8p8_Q2_2_2.5", {2, 2.5}}};
        xrangelims = {{"_x_0_0.25", {0.00001, 0.25}}};
        zrangelims = {{"_z_0.3_0.4", {0.3, 0.4}}, {"_z_0.4_0.60", {0.40, 0.60}}};
      }
      else if(tab=="numbers_x2_pim_8p8_hQ2_He3_kaon.txt" or tab=="numbers_x2_pip_8p8_hQ2_He3_kaon.txt"){
        Q2rangelims = {{"_8p8_Q2_3_8", {3, 8}}};
        xrangelims = {{"_x_0.25_0.7", {0.25,07}}};
        zrangelims = {{"_z_0.3_0.4", {0.3, 0.4}}, {"_z_0.4_0.60", {0.40, 0.60}}};
      }
      //std::cout<<" file name "<<tab<<std::endl;
      // Create directory
      std::string opath = ProcessedDataPath + "/" + ofolder;
      mkdir(opath.c_str(), ACCESSPERMS);
      
      // Reading table with YAML
      const YAML::Node exp = YAML::LoadFile(RawDataFolder + "/" + tab);

      // Create an input filestream
      std::ifstream rawfile(RawDataFolder + tab);

      // Make sure the file is open
      if(!rawfile.is_open()) throw std::runtime_error("Could not open file");

      std::string line;
      // Vector to hold the rows of data
      std::vector<DataRow> rows;
      if(rawfile.good())
      {
        while (std::getline(rawfile,line)){
          // If the line starts whith i, skip it.
          // (They are comments in the header of the rawdata file).
          if (line[0] == 'i')
            continue;
          rows.push_back(parseLine(line));
        }

      }
      

      // Conditions to separate the values in the columns into different files.
      // Loop on Q bin boundaries.
      for (auto const& Q2 : Q2rangelims)
      {
        //std::cout<<Q2.first<<"check"<<std::endl;
        for (auto const& x : xrangelims)
        {
          for (auto const& z : zrangelims)
          {
            //std::cout<<" check numbers x "<<x.second.first<<"  "<< x.second.second <<" z "<<z.second.first <<"  "<< z.second.second <<" Q2 "<< Q2.second.first <<"   "<< Q2.second.second<<std::endl;
            int charge =0;
            std::string name;
            if (rows[0].had=="K+") {charge=1;name="Kp";}
            else if (rows[0].had=="K-") {charge=-1;name="Km";}
            
            // File name cks up the m (i.e. Q) boundaries
            std::string ofileQ = "SoLID_Unpol" + Q2.first + x.first + z.first+"_"+name;

            // Initialize indexes vector for future selection
            // Vector to hold the indices of the selected rows
            std::vector<int> indexes;
            // Filter the rows based on the ranges
            for (const auto& row : rows) {
              //std::cout<<" check numbers"<<x.second.first<<" "<< row.x <<" "<< x.second.second <<" "<<z.second.first <<" "<<  row.z <<" "<< z.second.second <<" "<< Q2.second.first <<" "<< row.Q2 <<" "<< row.Q2 <<" "<< Q2.second.second<<std::endl;
              if (x.second.first < row.x && row.x < x.second.second && z.second.first < row.z && row.z < z.second.second && Q2.second.first < row.Q2 && row.Q2 < Q2.second.second) {
                indexes.push_back(row.index);
                //std::cout<<"check index"<<row.index<<std::endl;
              }
            }


            // Initialize result maps
            std::map<int, double> fdcross, fdstat, fdpT,fx,fy,fz,fQ2;

            // Prepare (outer) map for each output data file
            std::map<std::string, std::map<int, double>> filedata;

            // Finally select cross sections and their uncertainties.
            for (int i : indexes)
            {
              fdcross[i] = rows[i].value;
              fdstat[i]  = rows[i].stat;
              fdpT[i]    = rows[i].pT;
              fx[i]      = rows[i].x;
              fy[i]      = rows[i].y;
              fz[i]      = rows[i].z;
              fQ2[i]     = rows[i].Q2;
              //std::cout<<"check i "<<i<<" cross "<<fdcross[i]<<"pT"<<fdpT[i]<<std::endl;
            }

            filedata["cross"] = fdcross;
            filedata["stat"]  = fdstat;
            filedata["pT"]    = fdpT;
            filedata["x"]    = fx;
            filedata["y"]    = fy;
            filedata["z"]    = fz;
            filedata["Q2"]    = fQ2;

            // Plot labels
            std::map<std::string, std::string> labels
            {
              //{"xlabel", "#it{q}_{T} [GeV]"},
              //{"ylabel", "#frac{d^2#it{#sigma}}{d#it{Q}#d#it{q}_{T}} [pb GeV^{-2}]"},
              {"title", "SoLID, 3He target " + std::to_string(Q2.second.first) + " < Q2 < " + std::to_string(Q2.second.second)+","+ std::to_string(x.second.first) + " < x < " + std::to_string(x.second.second)+","+std::to_string(z.second.first) + " < z < " + std::to_string(z.second.second)},
                {"xlabelpy", "$P_{hT} \\rm{[GeV]}$"},
                {"ylabelpy", "$F^{K^"+std::string(1,rows[0].had[1])+"}_{\\rm UU,T }\\left(x, z, q_T, Q^2 \\right)$"},
                //{"xlabelpy", "$q_T \\rm{[GeV]}$"},
                //{"ylabelpy", "$\\frac{d^2\\sigma}{dQ dq_{T}}[\\rm{pb GeV}^{-2}]$"},
                {"titlepy", "SoLID unpol, " + std::to_string(Q2.second.first) + " < Q2 < " + std::to_string(Q2.second.second)}
            };

            /*
               NOTE on the conversion factor for the cross section.
               The raw data have a cross section expressed in cm**2/GeV**2/nucleon, but we would like to have it in NB./NUCLEON/GEV**2,
               then we have to convert: 1barn = 10**{-28}m**2 = 10**{-24}cm**2. Therefore, 1 cm**2 = 10**{24}barn= 10**{33}nb.

               NOTE on the calculation of y_min and y_max: y = arcsinh(sqrt{s}*xF/(2Q)).
               The value of x_min = 0, then y_min=0 for all bin in Q
               x_max = 1, then y_max = arcsinh(sqrt(s)/(2 Q_min)) for a specific Qmin < Q < Qmax bin
               */

            // Allocate emitter
            YAML::Emitter emit;

            // Write kinematics on the YAML emitter
            emit.SetFloatPrecision(8);
            emit.SetDoublePrecision(8);
            emit << YAML::BeginMap;
            emit << YAML::Key << "dependent_variables";
            emit << YAML::BeginSeq;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Value << YAML::Flow << labels;
            emit << YAML::Key << "qualifiers" << YAML::Value;
            emit << YAML::BeginSeq;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "process" << YAML::Key << "value" << YAML::Value << "SIDIS" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "observable" << YAML::Key << "value" << YAML::Value << "FUUT" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "target_isoscalarity" << YAML::Key << "value" << YAML::Value << "2/3" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "hadron" << YAML::Key << "value" << YAML::Value << "" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "charge" << YAML::Key << "value" << YAML::Value << charge << YAML::EndMap;
            //emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "beam" << YAML::Key << "value" << YAML::Value << rows[0].Eb << YAML::EndMap;
            //emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "hadron" << YAML::Key << "value" << YAML::Value << rows[0].had << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "prefactor" << YAML::Key << "value" << YAML::Value << 1 << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "Vs" << YAML::Key << "value" << YAML::Value << Vs << YAML::EndMap;
            //emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "Q2" << YAML::Key << "low" << YAML::Value << Q2.second.first << YAML::Key << "high" << YAML::Value << Q2.second.second  << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "Q" << YAML::Key << "low" << YAML::Value << sqrt(Q2.second.first) << YAML::Key << "high" << YAML::Value << sqrt(Q2.second.second)  << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            //emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "y" << YAML::Key
            //  << "low" << YAML::Value << asinh(0) << YAML::Key << "high" << YAML::Value << asinh(Vs / (2 * Q2.second.first)) << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "x" << YAML::Key
              << "low" << YAML::Value << x.second.first << YAML::Key << "high" << YAML::Value << x.second.second  << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "z" << YAML::Key
              << "low" << YAML::Value << z.second.first << YAML::Key << "high" << YAML::Value << z.second.second  << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "PS_reduction" << YAML::Key
               << "W" << YAML::Value << "2.3" << YAML::Key << "ymin" << YAML::Value << "0.01" << YAML::Key << "ymax" << YAML::Value << "0.95" << YAML::EndMap; 
               //<< "pTmin" << YAML::Value << "###" << YAML::Key << "etamin" << YAML::Value << "###" << YAML::Key << "etamax" << YAML::Value << "###" << YAML::EndMap; 
            emit << YAML::EndSeq;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& m : filedata["cross"])
            {
              //std::cout<<"check m"<<m.first<<std::endl;
              emit << YAML::BeginMap << YAML::Key << "errors" << YAML::Value << YAML::BeginSeq;

              // Oss: the conversion from cm**2 to barn is: 1b = 10**{-24}cm**2,
              //      but we need the cross section in pb: 1pb = 10**{-36}cm**2,
              //      <=> 1 cm**2 = 10**{36}pb
              // Remember conversion factor for the cross section
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "label" << YAML::Value << "unc" << YAML::Key << "value" << YAML::Value << (filedata["stat"][m.first])  << YAML::EndMap;
              //std::cout<<"check" <<filedata["stat"][m.first]<<std::endl;

              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "label" << YAML::Value << "add" << YAML::Key << "value" << YAML::Value << 0.1 << YAML::EndMap;
              emit << YAML::EndSeq;
              emit << YAML::Key << "value" << YAML::Value << (m.second) ;
              // emit << YAML::Key << "id"    << YAML::Value << m.first;
              emit << YAML::EndMap;
            }
            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::Key << "independent_variables";
            emit << YAML::BeginSeq;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Value << YAML::Flow;
            emit << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "PT" << YAML::Key << "units" << YAML::Value << "GEV" << YAML::EndMap;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            // for (auto const& p : filedata["pT"])
            //   {
            //     emit << YAML::Flow << YAML::BeginMap << YAML::Key << "value" << YAML::Value << p.second << YAML::Key << "high" << YAML::Value << (p.second + 0.125) << YAML::Key << "low" << YAML::Value << (p.second - 0.125) << YAML::Key
            //     /* Since Nangaparbat cross section is differential in qT, while this experimental set is differential in qT**2 and Q, we need to correct the value of the predictions by a factor DqT/(DQ * D(qT**2))*/
            //     << "factor" << YAML::Value << (p.second - p.first) / ((Q2.second.second - Q2.second.first) * (pow(p.second, 2) - pow(p.first, 2)))<< YAML::EndMap;
            //     // emit << YAML::Flow << YAML::BeginMap << YAML::Key << "value" << YAML::Value << p.second << YAML::EndMap;
            //   }
            for (auto const& p : filedata["pT"])
            {
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "value" << YAML::Value << p.second << YAML::Key
                /* Since Nangaparbat cross section is differential in qT, while this experimental set is differential in qT**2 and Q, we need to correct the value of the predictions by a factor DqT/(DQ * D(qT**2))*/
                << "factor" << YAML::Value << 1.0 / ((2*p.second)*(sqrt((Q2.second.second - Q2.second.first)/2)/((x.second.second-x.second.first)/2)/rows[0].Eb/0.93827208816))<< YAML::EndMap;
              //std::cout<<"check p.second"<<p.second<<std::endl;
              // emit << YAML::Flow << YAML::BeginMap << YAML::Key << "value" << YAML::Value << p.second << YAML::EndMap;
            }

            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "x" << YAML::EndMap;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& ix : filedata["x"])
              //std::cout<<" check x "<<ix.second<<std::endl;
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "high" << YAML::Value << x.second.second << YAML::Key << "low" << YAML::Value << x.second.first << YAML::Key << "value" << YAML::Value << ix.second << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "z" << YAML::EndMap;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& iz : filedata["z"])
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "high" << YAML::Value << z.second.second << YAML::Key << "low" << YAML::Value << z.second.first << YAML::Key << "value" << YAML::Value << iz.second << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "y" << YAML::EndMap;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& iy : filedata["y"])
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "high" << YAML::Value << 0.95 << YAML::Key << "low" << YAML::Value << 0.01 << YAML::Key << "value" << YAML::Value << iy.second << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::BeginMap;
            emit << YAML::Key << "header" << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "Q2" << YAML::EndMap;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& iQ2 : filedata["Q2"])
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "high" << YAML::Value << Q2.second.second << YAML::Key << "low" << YAML::Value << Q2.second.first << YAML::Key << "value" << YAML::Value << iQ2.second << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::EndMap;
            emit << YAML::EndSeq;
            emit << YAML::EndMap;

            // Dump table to file
            std::ofstream fout(opath + "/" + ofileQ + ".yaml");
            fout << emit.c_str() << std::endl;
            fout.close();

            filenames.push_back(ofileQ);
          }

        }//loop over Q2
      }//loop overr x
    }//loop over z
       
    // Map to space properly the dataset names in datasets.yaml
    std::map<int, std::string> spaces = {{51," "}, {50,"  "}, {49,"   "}, {48, "    "}, {47, "     "}, {46, "      "},{45, "       "},{44,"        "}};

    // Produce outputnames to put in datasets.yaml
    std::string outputnames;
    for (std::string name : filenames)
      outputnames += "  - {name: " + name + "," + spaces[name.size()] + "file: " + name + ".yaml}\n";

    return outputnames;


    //return
    //  "check";
  }
}

int main() {
  std::string raw_data_path = ".";
  std::string processed_data_path = ".";
  //bool pdf_error = false;
  // Dataset file created dunring the filtering
  std::ofstream fout(processed_data_path + "/datasets.yaml");

  fout << "#SoLID 11GeV Pion:\n";

  fout<< NangaParbat::PreprocessSoLID(raw_data_path, processed_data_path);
  
  fout << "#SoLID 8.8GeV Pion:\n";

  fout<< NangaParbat::PreprocessSoLID8p8(raw_data_path, processed_data_path);
  
  fout << "#SoLID 11GeV Kaon:\n";

  fout<< NangaParbat::PreprocessSoLIDKaon(raw_data_path, processed_data_path);
  
  fout << "#SoLID 8.8GeV Kaon:\n";

  fout<< NangaParbat::PreprocessSoLID8p8Kaon(raw_data_path, processed_data_path);

  return 0;
}

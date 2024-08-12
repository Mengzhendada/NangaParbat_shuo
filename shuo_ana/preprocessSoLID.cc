/*
 * Authors: Simone Venturini
 *          Chiara Bissolotti: chiara.bissolotti01@universitadipavia.it
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
    row.x = std::stod(item);


    // Skip the next 2 columns (xl,xr)
    std::getline(ss, item, ','); // Skip xr
    std::getline(ss, item, ','); // Skip xl

    // Read y
    std::getline(ss, item, ',');
    row.y = std::stod(item);

    // Read z
    std::getline(ss, item, ',');
    row.z = std::stod(item);

    // Skip the next 2 columns (zl, zr)
    std::getline(ss, item, ','); // Skip zl
    std::getline(ss, item, ','); // Skip zr

    // Read Q2
    std::getline(ss, item, ',');
    row.Q2 = std::stod(item);

    // Skip the next 2 columns (Q2l, Q2r)
    std::getline(ss, item, ','); // Skip Q2l
    std::getline(ss, item, ','); // Skip Q2r

    // Read pT
    std::getline(ss, item, ','); 
    row.pT = std::stod(item);

    // Skip the next 3 columns (pTl, pTr, obs)
    std::getline(ss, item, ','); // Skip pTl
    std::getline(ss, item, ','); // Skip pTr
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
  std::string PreprocessSoLID(std::string const& RawDataPath, std::string const& ProcessedDataPath, bool const& PDFError)
  {
    std::cout << "Processing SoLID data ..." << std::endl;

    // Path to the raw-data folder
    const std::string RawDataFolder = RawDataPath + "/SoLID_Unpol/";


    // Vector of tables to process
    const std::vector<std::string> tables = {"numbers_x1_pip_11.txt"};

    // Output folder
    const std::string ofolder = "oSoLID_Unpol";

    // Vector for output names
    std::vector<std::string> filenames;

    // Initialize naming map for the Q-integration ranges (the first element is for the name of the output data file)
    const std::map<std::string, std::pair<double, double>> Q2rangelims = {{"_Q2_1_1.50", {1.0, 1.50}}, {"_Q2_1.50_2.0", {1.50, 2.0}}, {"_Q2_2.0_2.5", {2.0, 2.5}}};
    const std::map<std::string, std::pair<double, double>> xrangelims = {{"_x_0_0.25", {0, 0.25}}, {"_x_0.25_0.50", {0.250, 0.50}}};
    const std::map<std::string, std::pair<double, double>> zrangelims = {{"_z_0.3_0.35", {0.3, 0.35}}, {"_z_0.35_0.40", {0.350, 0.40}},{"_z_0.40_0.45",{0.40,0.45}},{"_z_0.45_0.5", {0.45, 0.5}}, {"_z_0.5_0.55", {0.50, 0.55}},{"_z_0.55_0.6",{0.55,0.6}}};

    // Vs ( = sqrt(2*M*E) for fixed target experiments, here E = 11 GeV)
    const double Vs = 4.7;

    // Create directory
    std::string opath = ProcessedDataPath + "/" + ofolder;
    mkdir(opath.c_str(), ACCESSPERMS);

    // Loop over tables
    for (auto const& tab : tables)
    {
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
          // If the line starts whith #, skip it.
          // (They are comments in the header of the rawdata file).
          if (line[0] == '#')
            continue;
          rows.push_back(parseLine(line));
        }

      }
      

      // Conditions to separate the values in the columns into different files.
      // Loop on Q bin boundaries.
      for (auto const& Q2 : Q2rangelims)
      {
        std::cout<<Q2.first<<"check"<<std::endl;
        for (auto const& x : xrangelims)
        {
          for (auto const& z : zrangelims)
          {
            // File name picks up the m (i.e. Q) boundaries
            std::string ofileQ = "SoLID_Unpol" + Q2.first + x.first + z.first;

            // Initialize indexes vector for future selection
            // Vector to hold the indices of the selected rows
            std::vector<int> indexes;
            // Filter the rows based on the ranges
            for (const auto& row : rows) {
              std::cout<<" check numbers"<<x.second.first<<" "<< row.x <<" "<< x.second.second <<" "<<z.second.first <<" "<<  row.z <<" "<< z.second.second <<" "<< Q2.second.first <<" "<< row.Q2 <<" "<< row.Q2 <<" "<< Q2.second.second<<std::endl;
              if (x.second.first < row.x && row.x < x.second.second && z.second.first < row.z && row.z < z.second.second && Q2.second.first < row.Q2 && row.Q2 < Q2.second.second) {
                indexes.push_back(row.index);
                std::cout<<"check index"<<row.index<<std::endl;
              }
            }


            // Initialize result maps
            std::map<int, double> fdcross, fdstat, fdpT;

            // Prepare (outer) map for each output data file
            std::map<std::string, std::map<int, double>> filedata;

            // Finally select cross sections and their uncertainties.
            for (int i : indexes)
            {
              fdcross[i] = rows[i].value;
              fdstat[i]  = rows[i].stat;
              fdpT[i]    = rows[i].pT;
              std::cout<<"i "<<i<<" cross "<<fdcross[i]<<"pT"<<fdpT[i]<<std::endl;
            }

            filedata["cross"] = fdcross;
            filedata["stat"]  = fdstat;
            filedata["pT"]    = fdpT;



            // Plot labels
            std::map<std::string, std::string> labels
            {
              //{"xlabel", "#it{q}_{T} [GeV]"},
              //{"ylabel", "#frac{d^2#it{#sigma}}{d#it{Q}#d#it{q}_{T}} [pb GeV^{-2}]"},
              {"title", "SoLID, 3He target " + std::to_string(Q2.second.first) + " < Q2 < " + std::to_string(Q2.second.second)+","+ std::to_string(x.second.first) + " < x < " + std::to_string(x.second.second)+","+std::to_string(z.second.first) + " < z < " + std::to_string(z.second.second)},
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
            //emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "target" << YAML::Key << "value" << YAML::Value << 0.4025 << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "beam" << YAML::Key << "value" << YAML::Value << rows[0].Eb << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "hadron" << YAML::Key << "value" << YAML::Value << rows[0].had << YAML::EndMap;
            //emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "prefactor" << YAML::Key << "value" << YAML::Value << 1 << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "Vs" << YAML::Key << "value" << YAML::Value << Vs << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "Q2" << YAML::Key
              << "low" << YAML::Value << Q2.second.first << YAML::Key << "high" << YAML::Value << Q2.second.second  << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            //emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "y" << YAML::Key
            //  << "low" << YAML::Value << asinh(0) << YAML::Key << "high" << YAML::Value << asinh(Vs / (2 * Q2.second.first)) << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "x" << YAML::Key
              << "low" << YAML::Value << x.second.first << YAML::Key << "high" << YAML::Value << x.second.second  << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "z" << YAML::Key
              << "low" << YAML::Value << z.second.first << YAML::Key << "high" << YAML::Value << z.second.second  << YAML::Key << "integrate" << YAML::Value << "true" << YAML::EndMap;
            /* emit << YAML::Flow << YAML::BeginMap << YAML::Key << "name" << YAML::Value << "PS_reduction" << YAML::Key
               << "pTmin" << YAML::Value << "###" << YAML::Key << "etamin" << YAML::Value << "###" << YAML::Key << "etamax" << YAML::Value << "###" << YAML::EndMap; */
            emit << YAML::EndSeq;
            emit << YAML::Key << "values" << YAML::Value;
            emit << YAML::BeginSeq;
            for (auto const& m : filedata["cross"])
            {
              std::cout<<"m"<<m.first<<std::endl;
              emit << YAML::BeginMap << YAML::Key << "errors" << YAML::Value << YAML::BeginSeq;

              // Oss: the conversion from cm**2 to barn is: 1b = 10**{-24}cm**2,
              //      but we need the cross section in pb: 1pb = 10**{-36}cm**2,
              //      <=> 1 cm**2 = 10**{36}pb
              // Remember conversion factor for the cross section
              emit << YAML::Flow << YAML::BeginMap << YAML::Key << "label" << YAML::Value << "unc" << YAML::Key << "value" << YAML::Value << (filedata["stat"][m.first])  << YAML::EndMap;
              std::cout<<"check" <<filedata["stat"][m.first]<<std::endl;

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
              std::cout<<"check p.second"<<p.second<<std::endl;
              // emit << YAML::Flow << YAML::BeginMap << YAML::Key << "value" << YAML::Value << p.second << YAML::EndMap;
            }


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
        /*
        // Map to space properly the dataset names in datasets.yaml
        std::map<int, std::string> spaces = {{18, " "}, {17, "  "}, {16, "   "}};

        // Produce outputnames to put in datasets.yaml
        std::string outputnames;
        for (std::string name : filenames)
        outputnames += "  - {name: " + name + "," + spaces[name.size()] + "file: " + name + ".yaml}\n";

        return outputnames;
        */
        return
          "check";
        // "  - {name: E615_Q_4.05_4.50,    file: E615_Q_4.05_4.50.yaml}\n"
        // "  - {name: E615_Q_4.50_4.95,    file: E615_Q_4.50_4.95.yaml}\n"
        // "  - {name: E615_Q_4.95_5.40,    file: E615_Q_4.95_5.40.yaml}\n"
        // "  - {name: E615_Q_5.40_5.85,    file: E615_Q_5.40_5.85.yaml}\n"
        // "  - {name: E615_Q_5.85_6.75,    file: E615_Q_5.85_6.75.yaml}\n"
        // "  - {name: E615_Q_6.75_7.65,    file: E615_Q_6.75_7.65.yaml}\n"
        // "  - {name: E615_Q_7.65_9.00,    file: E615_Q_7.65_9.00.yaml}\n"
        // "#  - {name: E615_Q_9.00_10.35,   file: E615_Q_9.00_10.35.yaml}\n"
        // "#  - {name: E615_Q_10.35_11.70,  file: E615_Q_10.35_11.70.yaml}\n"
        // "  - {name: E615_Q_11.70_13.05,  file: E615_Q_11.70_13.05.yaml}\n";
  }
}

int main() {
  std::string raw_data_path = ".";
  std::string processed_data_path = ".";
  bool pdf_error = false;

  std::string result = NangaParbat::PreprocessSoLID(raw_data_path, processed_data_path, pdf_error);
  std::cout << result << std::endl;

  return 0;
}

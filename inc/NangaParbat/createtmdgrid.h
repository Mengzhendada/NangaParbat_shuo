//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "NangaParbat/parameterisation.h"

#include <vector>
#include <memory>
#include <yaml-cpp/yaml.h>
#include <apfel/apfelxx.h>

namespace NangaParbat
{
  /**
   * @brief Structure that contains the 3D grid in Q, x, and qT.
   */
  struct ThreeDGrid
  {
    std::vector<double> Qg;       //!< Grid in Q
    std::vector<double> xg;       //!< Grid in x
    std::vector<double> qToQg;    //!< Grid in  qT/Q
  };

  /**
   * @brief Function that returns ThreeDGrid object better tuned for
   * PDFs or FFs, according to the input string.
   */
  ThreeDGrid InterGrid(std::string const& pf)
  {
    ThreeDGrid grid{};
    if (pf == "pdf")
      {
        grid.Qg =
        {
          1.000000e+00, 1.118034e+00, 1.224745e+00, 1.400000e+00, 1.581139e+00, 1.788854e+00,
          2.000000e+00, 2.236068e+00, 2.529822e+00, 2.828427e+00, 3.162278e+00, 3.464102e+00,
          4.750000e+00, 5.099020e+00, 6.324555e+00, 7.100000e+00, 8.000000e+00, 1.000000e+01,
          1.118034e+01, 1.224745e+01, 1.400000e+01, 1.581139e+01, 1.788854e+01, 2.000000e+01,
          2.236068e+01, 2.529822e+01, 2.828427e+01, 3.162278e+01, 3.464102e+01, 4.750000e+01,
          5.099020e+01, 6.324555e+01, 7.100000e+01, 8.000000e+01, 1.000000e+02, 1.118034e+02,
          1.224745e+02, 1.400000e+02, 1.581139e+02, 1.788854e+02, 2.000000e+02
        };
        grid.xg =
        {
          1.000000e-05, 2.000000e-05, 4.000000e-05, 6.000000e-05, 8.000000e-05,
          1.000000e-04, 2.000000e-04, 4.000000e-04, 6.000000e-04, 8.000000e-04,
          1.000000e-03, 1.500000e-03, 2.000000e-03, 2.500000e-03, 3.000000e-03,
          3.500000e-03, 4.000000e-03, 4.500000e-03, 5.000000e-03, 5.500000e-03,
          6.000000e-03, 6.500000e-03, 7.000000e-03, 7.500000e-03, 8.000000e-03,
          8.500000e-03, 9.000000e-03, 9.250000e-03, 9.500000e-03, 9.750000e-03,
          1.000000e-02, 1.500000e-02, 2.000000e-02, 2.500000e-02, 3.000000e-02,
          3.500000e-02, 4.000000e-02, 4.500000e-02, 5.000000e-02, 5.500000e-02,
          6.000000e-02, 6.500000e-02, 7.000000e-02, 7.500000e-02, 8.000000e-02,
          8.500000e-02, 9.000000e-02, 9.250000e-02, 9.500000e-02, 9.750000e-02,
          1.000000e-01, 1.500000e-01, 2.000000e-01, 2.500000e-01, 3.000000e-01,
          3.500000e-01, 4.000000e-01, 4.500000e-01, 5.000000e-01, 5.500000e-01,
          6.000000e-01, 6.500000e-01, 7.000000e-01, 7.500000e-01, 8.000000e-01,
          8.500000e-01, 9.000000e-01, 9.250000e-01, 9.500000e-01, 9.750000e-01,
          1.000000e+00
        };
        grid.qToQg =
        {
          0.0001, 0.0010, 0.0025, 0.0050, 0.0075, 0.0100, 0.0200, 0.0300, 0.0400,
          0.0500, 0.0600, 0.0700, 0.0800, 0.0900, 0.1000, 0.1250, 0.1500, 0.1750,
          0.2000, 0.2250, 0.2500, 0.2750, 0.3000, 0.3500, 0.4000, 0.4500, 0.5000,
          0.5500, 0.6000, 0.6500, 0.7000, 0.8000, 0.9000, 1
        };
      }
    else
      {
        grid.Qg =
        {
          1.000000e+00, 1.118034e+00, 1.224745e+00, 1.400000e+00, 1.581139e+00, 1.788854e+00,
          2.000000e+00, 2.236068e+00, 2.529822e+00, 2.828427e+00, 3.162278e+00, 3.464102e+00,
          4.750000e+00, 5.099020e+00, 6.324555e+00, 7.100000e+00, 8.000000e+00, 1.000000e+01,
          1.118034e+01, 1.224745e+01, 1.400000e+01, 1.581139e+01, 1.788854e+01, 2.000000e+01,
          2.236068e+01, 2.529822e+01, 2.828427e+01, 3.162278e+01, 3.464102e+01, 4.750000e+01,
          5.099020e+01, 6.324555e+01, 7.100000e+01, 8.000000e+01, 1.000000e+02, 1.118034e+02,
          1.224745e+02, 1.400000e+02, 1.581139e+02, 1.788854e+02, 2.000000e+02
        };
        grid.xg = std::vector<double>
        {
          1.000000e-02, 1.200000e-02, 1.500000e-02, 1.700000e-02, 2.000000e-02, 2.200000e-02, 2.500000e-02, 3.000000e-02,
          3.500000e-02, 3.700000e-02, 4.000000e-02, 4.200000e-02, 4.500000e-02, 4.700000e-02, 5.000000e-02, 5.500000e-02,
          6.000000e-02, 6.200000e-02, 6.500000e-02, 6.700000e-02, 7.000000e-02, 7.200000e-02, 7.500000e-02, 7.700000e-02,
          8.000000e-02, 8.200000e-02, 8.500000e-02, 8.700000e-02, 9.000000e-02, 9.200000e-02, 9.500000e-02, 9.700000e-02,
          1.000000e-01, 1.200000e-01, 1.500000e-01, 1.700000e-01, 2.000000e-01, 2.200000e-01, 2.500000e-01, 3.000000e-01,
          3.500000e-01, 3.700000e-01, 4.000000e-01, 4.200000e-01, 4.500000e-01, 4.700000e-01, 5.000000e-01, 5.500000e-01,
          6.000000e-01, 6.200000e-01, 6.500000e-01, 6.700000e-01, 7.000000e-01, 7.200000e-01, 7.500000e-01, 7.700000e-01,
          8.000000e-01, 8.200000e-01, 8.500000e-01, 8.700000e-01, 9.000000e-01, 9.200000e-01, 9.500000e-01, 9.700000e-01,
          1.000000e+00
        };
        grid.qToQg = std::vector<double>
        {
          0.0001,
          0.0010, 0.0025, 0.0050, 0.0075, 0.0100, 0.0200, 0.0300, 0.0400,
          0.0500, 0.0600, 0.0700, 0.0800, 0.0900, 0.1000, 0.1100, 0.1200,
          0.1300, 0.1400, 0.1500, 0.1600, 0.1700, 0.1800, 0.1900, 0.2000,
          0.2200, 0.2400, 0.2600, 0.2800, 0.3000, 0.3200, 0.3400, 0.3600,
          0.3800, 0.4000, 0.4500, 0.5000, 0.5500, 0.6000, 0.6500, 0.7000,
          0.8000, 0.9000, 1
        };
      }
    return grid;
  };

  /**
   * @brief Map of perturbative orders
   */
  const std::map<int, std::string> PtOrderMap{{0, "LL"}, {1, "NLL"}, {-1, "NLL'"}, {2, "NNLL"}, {-2, "NNLL'"}, {3, "NNNLL"}};

  /**
   * @brief This function encapsulates and streamlines the production
   * of an interpolation grid strting from the report produced by a
   * NangaParbat fit.
   * @param ReportFolder: path to the report folder
   * @param Output: name of the output grid
   * @param pf: whether PDFs ("pdf") of FFs ("ff")
   */
  void ProduceTMDGrid(std::string const& ReportFolder, std::string const& Output, std::string const& distype = "pdf");

  /**
   * @brief Function that produces the TMD interpolation grid in
   * momentum space. This is supposed to resamble an LHAPDF grid for
   * the TMDs. We use plain YAML format.
   * @param config: the YAML node with the theory settings
   * @param parameterisation: the parameterisation type
   * @param params: the vector of parameters to be used for the tabulation
   * @param pf: whether PDFs ("pdf") of FFs ("ff")
   * @return a YAML emitter
   */
  std::unique_ptr<YAML::Emitter> EmitTMDGrid(YAML::Node          const& config,
                                             std::string         const& parameterisation,
                                             std::vector<double> const& params,
                                             std::string         const& pf,
                                             ThreeDGrid          const& tdg);

  /**
   * @brief Function that produces the info file of the TMD set. This
   * is suppose to resamble an LHAPDF info file for the TMDs. We use
   * plain YAML format.
   * @param config: the YAML not with the theory settings
   * @param pf: whether PDFs ("pdf") of FFs ("ff")
   * @return a YAML emitter
   */
  std::unique_ptr<YAML::Emitter> EmitTMDInfo(YAML::Node  const& config,
                                             int         const& NumMembers,
                                             std::string const& pf,
                                             std::string const& SetDesc = "Set produced with NangaParbat + APFEL++",
                                             std::string const& Authors = "A. Bacchetta, V. Bertone, C. Bissolotti, G. Bozzi, F. Delcarro, F. Piacenza, M. Radici",
                                             std::string const& Reference = "arXiv:xxxx.xxxxx",
                                             std::string const& SetIndex = "000000",
                                             std::string const& Format = "TMDlib1",
                                             std::string const& DataVersion = "1",
                                             std::string const& ErrorType = "Monte Carlo");
}
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "NangaParbat/cutfactory.h"
#include "NangaParbat/zcut.h"
#include "NangaParbat/xcut.h"
#include "NangaParbat/qcut.h"
#include "NangaParbat/xzcut.h"
#include "NangaParbat/wcut.h"

namespace NangaParbat
{
  //_________________________________________________________________________________
  std::shared_ptr<Cut> CutFactory::GetInstance(DataHandler const& dataset, std::string const& name, double const& min, double const& max, std::vector<double> const& params)
  {
    if (name == "zcut")
      return std::shared_ptr<NangaParbat::Cut> {new NangaParbat::ZCut{dataset, min, max}};
    else if (name == "xcut")
      return std::shared_ptr<NangaParbat::Cut> {new NangaParbat::XCut{dataset, min, max}};
    else if (name == "Qcut")
      return std::shared_ptr<NangaParbat::Cut> {new NangaParbat::QCut{dataset, min, max}};
    else if (name == "xzcut")
      return std::shared_ptr<NangaParbat::Cut> {new NangaParbat::XZCut{dataset, min, max, params[0]}};
    else if (name == "xzcutQuad")
      return std::shared_ptr<NangaParbat::Cut> {new NangaParbat::XZCutQuadratic{dataset, min, max}};
    else if (name == "xzcutLog")
      return std::shared_ptr<NangaParbat::Cut> {new NangaParbat::XZCutLog{dataset, min, max, params[0]}};
    else if (name == "Wcut")
      return std::shared_ptr<NangaParbat::Cut> {new NangaParbat::WCut{dataset, min, max}};
    else
      throw std::runtime_error("[CutFactory::GetInstance]: Unknown cut.");
  }
}

//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Alessandro Bacchetta: alessandro.bacchetta@unipv.it

#pragma once

#include "NangaParbat/parameterisation.h"

#include <math.h>

namespace NangaParbat
{
  /**
   * @brief Parameterisation derived from the "Parameterisation"
   * mother class used to fit the pion TMDs.
   */
  class MAPTMDPion3: public NangaParbat::Parameterisation
  {
  public:

    MAPTMDPion3():
      Parameterisation{"MAPTMDPion3", 2, std::vector<double> {0.2343, 0.2788, 0.8466, 1.3806, 1.9872, 0.1590, 0.4534, 3.843, 0.0276, 0.00392, 12.551, 4.141, 0.2343, 0.2788 ,0.6089}} { };

    double Evaluate(double const& x, double const& b, double const& zeta, int const& ifunc) const
    {
      if (ifunc < 0 || ifunc >= this->_nfuncs)
        throw std::runtime_error("[MAPTMDPion3::Evaluate]: function index out of range");

      // If the value of 'x' exceeds one returns zero
      if (x >= 1)
        return 0;

      // Evolution
      const double g2   = this->_pars[0];
      const double evol = exp( - pow(g2, 2) * log(zeta / _Q02) * b * b / 4 );

      if (ifunc == 0)   // TMD PDFs proton MAP22
        {
          const double N1      = this->_pars[1];
          const double sigma1  = this->_pars[2];
          const double alpha1  = this->_pars[3];
          const double lambda  = this->_pars[4];
          const double N1B     = this->_pars[5];
          const double sigma2  = this->_pars[6];
          const double alpha2  = this->_pars[7];
          const double lambda2 = this->_pars[8];
          const double N1C     = this->_pars[9];
          const double sigma3  = this->_pars[10];
          const double alpha3  = this->_pars[11];
          const double xhat    = 0.1;
          const double g1      = N1 * pow(x / xhat, sigma1) * pow((1 - x) / (1 - xhat), pow(alpha1, 2));
          const double g1B     = N1B * pow(x / xhat, sigma2) * pow((1 - x) / (1 - xhat), pow(alpha2, 2));
          const double g1C     = N1C * pow(x / xhat, sigma3) * pow((1 - x) / (1 - xhat), pow(alpha3, 2));
          return evol * ( g1 * exp( - g1 * pow(b / 2, 2))
                          +  pow(lambda, 2)  * pow(g1B, 2) * ( 1 - g1B * pow(b / 2, 2)) * exp( - g1B * pow(b / 2, 2)) + g1C * pow(lambda2, 2) * exp( - g1C * pow(b / 2, 2)))
                 / ( g1 +  pow(lambda, 2)  * pow(g1B, 2) + g1C * pow(lambda2, 2));
        }
      else    // TMD PDFs pion
        {
          const double N1pi     = this->_pars[12];
          const double sigmapi  = this->_pars[13];
          const double alphapi  = this->_pars[14];
          const double xhat    = 0.1;
          const double g1pi      = N1pi * pow(x / xhat, sigmapi) * pow((1 - x) / (1 - xhat), pow(alphapi, 2));
          return evol  * exp( - g1pi * pow( b / 2, 2) );
        }

    };

    std::string LatexFormula() const
    {
      std::string formula;
      formula  = R"delimiter($$f_{\rm NP}(x,\zeta, b_T)= \exp \left( g_{1\pi}(x) b_T^2 / 4 \right) \exp\left[- g^2_2 \log\left(\frac{\zeta}{Q_0^2}\right) b_T^2/4 \right]$$)delimiter"; 
      formula += R"delimiter($$g_{1\pi}(x) = N_{1\pi} \frac{x^{\sigma_{\pi}}(1-x)^{\alpha^2_{\pi}}}{\hat{x}^{\sigma_{\pi}}(1-\hat{x})^{\alpha^2_{\pi}}}$$)delimiter";
      formula += R"delimiter($$Q_0^2 = 1\;{\rm GeV}^2$$)delimiter";
      formula += R"delimiter($$\hat{x} = 0.1$$)delimiter";
      return formula;
    };

    std::vector<std::string> GetParameterNames() const
    {
      return {R"delimiter($g_2$)delimiter",
              R"delimiter($N_1$)delimiter",
              R"delimiter($\sigma1$)delimiter",
              R"delimiter($\alpha1$)delimiter",
              R"delimiter($\lambda$)delimiter",
              R"delimiter($N_1B$)delimiter",
              R"delimiter($\sigma2$)delimiter",
              R"delimiter($\alpha2$)delimiter",
              R"delimiter($\lambda2$)delimiter",
              R"delimiter($N1C$)delimiter",
              R"delimiter($\sigma3$)delimiter",
              R"delimiter($\alpha3$)delimiter",
              R"delimiter($N_{1\pi}$)delimiter",
              R"delimiter($\sigma_\pi$)delimiter",
              R"delimiter($\alpha_\pi$)delimiter"};
    };

    std::string GetDescription() const
    {
      return "Parameterisation used for the Pavia 2022 TMD analysis.";
    };

  private:
    double              const _Q02 = 1;
  };
}

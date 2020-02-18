// From the STL:
#include <iostream>
#include <iomanip>
#include <limits>

using namespace std;

// From bpp-core:
#include <Bpp/Version.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RASTools.h>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/MixedSubstitutionModel.h>
#include <Bpp/Phyl/Model/Protein/CoalaCore.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/FrequenciesSet/MvaFrequenciesSet.h>
#include <Bpp/Phyl/Model/FrequenciesSet/FrequenciesSet.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/BppOFrequenciesSetFormat.h>


using namespace bpp;
using namespace std;

#include <map>
#include <vector>

/******************************************************************************/
/**************************** Auxiliary functions *****************************/
/******************************************************************************/


void printModelParameters(DiscreteRatesAcrossSitesTreeLikelihood* tl)
  {
    double logL = tl->getValue();
    ApplicationTools::displayResult("Log likelihood", TextTools::toString(-logL, 15));
    ParameterList parameters = tl->getSubstitutionModelParameters();
    for (size_t i = 0; i < parameters.size(); i++)
    {
      ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
    }
    parameters = tl->getRateDistributionParameters();
    for (size_t i = 0; i < parameters.size(); i++)
    {
      ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
    }
    ParameterList pl = tl->getBranchLengthsParameters();
    for(unsigned int i = 0; i < pl.size(); i++)
    {
       ApplicationTools::displayResult(pl[i].getName(), TextTools::toString(pl[i].getValue())); 
    }
}

VectorSiteContainer* process_alignment(Alphabet* alphabet, BppApplication bppml)
{
    VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, bppml.getParams()); // here, gaps will be converted to unknown characters
    VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, bppml.getParams(), "", true, false); // convert the alignemnt to condensed format of unique sites
    delete allSites; // delete the non-condenced intance of the sequence data
    return(sites);

}

Tree* process_tree(BppApplication bppml)
{
    Tree* tree = PhylogeneticsApplicationTools::getTree(bppml.getParams());
    TreeTemplate<Node> ttree(*tree);
    vector<Node *> nodes = ttree.getNodes();
    for(unsigned int i = 0; i < nodes.size(); i++)
    {
      if(nodes[i]->isLeaf())
        nodes[i]->setName("Leaf" + TextTools::toString(nodes[i]->getId()));
      else
        nodes[i]->setName("Internal" + TextTools::toString(nodes[i]->getId()));
      ApplicationTools::displayResult(nodes[i]->getName(), TextTools::toString(nodes[i]->getId())); 
    }
    // try to add father to Internal4 (parent of Leaf3) - no such function - need to use detailsSimulation instance

    string initBrLenMethod = ApplicationTools::getStringParameter("init.brlen.method", bppml.getParams(), "Input", "", true, 1); // process tree branches lengths
    string cmdName;
    map<string, string> cmdArgs;
    KeyvalTools::parseProcedure(initBrLenMethod, cmdName, cmdArgs); // this line process cmdName - which dictates weather the tree should be processed as is, or ultrameterized
    if (cmdName == "Clock") // if needed, turn the tree into an ultrametric tree
    {
      TreeTools::convertToClockTree(*tree, tree->getRootId(), true);
    }
    return(tree);
}

/******************************************************************************/
/*********************************** Main *************************************/
/******************************************************************************/

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*       Test non-homogenous models    " << BPP_VERSION << "      *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  try
  { 

    // process input from params file
    BppApplication bppml(args, argv, "bppML");

    map<string,string> parans = bppml.getParams(); // debug

    /* process alphabet type */
    Alphabet* alphabet = SequenceApplicationTools::getAlphabet(bppml.getParams(), "", false);
    unique_ptr<GeneticCode> gCode;
    CodonAlphabet* codonAlphabet = dynamic_cast<CodonAlphabet*>(alphabet);
    if (codonAlphabet) {
      string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppml.getParams(), "Standard", "", true, true);
      ApplicationTools::displayResult("Genetic Code", codeDesc);
      gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
    }


    /* process alignment */
    VectorSiteContainer* sites = process_alignment(alphabet, bppml);

    /* process tree */
    Tree* tree = process_tree(bppml);
    
    /* test a double non homogenous model */
    cout << "\nDefine branch model - partition 1 " << endl;
    bppml.getParam("nonhomogeneous")="general";
    bppml.getParam("nonhomogeneous.number_of_models") = "2";
    // define the exact same model twice - and alias that parameters
    bppml.getParam("model1") ="YNGP_M2(kappa=1,omega_0=0.5,omega_2=2.0,frequencies=F3X4,theta1=0.333333,theta2=0.333333,initFreqs=observed,initFreqs.observedPseudoCount=1)";
    bppml.getParam("model2") = "YNGP_M2(kappa=YNGP_M2.kappa_1,omega_0=YNGP_M2.omega_0_1,omega_2=YNGP_M2.omega_2_1,theta1=YNGP_M2.theta1_1,theta2=YNGP_M2.theta2_1,frequencies=F3X4,initFreqs=observed,initFreqs.observedPseudoCount=1)";
    bppml.getParam("nonhomogeneous.stationarity") = "yes"; // constrain root frequencies to be the same as stationary (since RELAX is a time reversible model, this should not cause issues)


    bppml.getParam("site.number_of_paths") = "2";
    bppml.getParam("site.path1") = "model1[YNGP_M2.omega_0]&model2[YNGP_M2.omega_0]";
    bppml.getParam("site.path2") = "model1[YNGP_M2.omega_2]&model2[YNGP_M2.omega_2]";

    // generate the branch-site model
    SubstitutionModelSet* modelSet_1 = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, gCode.get(), sites, bppml.getParams());
    DiscreteDistribution* bs_rDist_1 = PhylogeneticsApplicationTools::getRateDistribution(bppml.getParams());
    DiscreteRatesAcrossSitesTreeLikelihood* bs_tl_1 = new RNonHomogeneousMixedTreeLikelihood(*tree, *sites, dynamic_cast<MixedSubstitutionModelSet*>(modelSet_1), bs_rDist_1, true, true);
    bs_tl_1->initialize();

    /* compute likelihood */
    cout << "\nCompute log likelihood" << endl;
    // should output the same log likelihood: -2215.84197049057 instead of -2216.76196291407
    printModelParameters(bs_tl_1); // debug - print model parameters

    // clasiffy branches into BG (model1) and FB (model2)
    cout << "\nDefine branch model - partition 2 " << endl;
    bppml.getParam("model1.nodes_id") = "0,1,2,3";
    bppml.getParam("model2.nodes_id") = "4,5,6,7";
    // generate the branch-site model
    SubstitutionModelSet* modelSet_2 = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, gCode.get(), sites, bppml.getParams());
    DiscreteDistribution* bs_rDist_2 = PhylogeneticsApplicationTools::getRateDistribution(bppml.getParams());
    DiscreteRatesAcrossSitesTreeLikelihood* bs_tl_2 = new RNonHomogeneousMixedTreeLikelihood(*tree, *sites, dynamic_cast<MixedSubstitutionModelSet*>(modelSet_2), bs_rDist_2, true, true);
    bs_tl_2->initialize();

    /* compute likelihood */
    cout << "\nCompute log likelihood" << endl;
    // should output the same log likelihood: -2215.84197049057 instead of -2216.76196291407
    printModelParameters(bs_tl_2); // debug - print model parameters

    /* free resources */
    delete alphabet;
    delete sites;
    delete tree;
    if (modelSet_1) delete modelSet_1;
    delete bs_rDist_1;
    delete bs_tl_1;
    if (modelSet_2) delete modelSet_2;
    delete bs_rDist_2;
    delete bs_tl_2;
    bppml.done();

  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}
// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/* \file Zoltan2_MetricAnalyzer.hpp
 * \brief Used by the Zoltan2 test driver for running \
          simple pass fail tests based on ranges of problem metrics.
 */
#ifndef ZOLTAN2_METRIC_ANALYZER_HPP
#define ZOLTAN2_METRIC_ANALYZER_HPP

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_Typedefs.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Teuchos_FileInputSource.hpp>

#include <sstream>
#include <string>
#include <iostream>

using Teuchos::ParameterList;
using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using namespace Zoltan2_TestingFramework;

// these define key names which convert to an API call
#define API_STRING_getWeightImbalance "imbalance"
#define API_STRING_getTotalEdgeCuts "total edge cuts"
#define API_STRING_getMaxEdgeCuts "max edge cuts"

// these define the options for a particular type
#define KEYWORD_PARAMETER_NAME "check" // would usually be the first entry and identify the API call
#define WEIGHT_PARAMETER_NAME "weight"
#define UPPER_PARAMETER_NAME "upper"
#define LOWER_PARAMETER_NAME "lower"
#define NORMED_PARAMETER_NAME "normed"

#define UNDEFINED_PARAMETER_INT_INDEX -1	// didn't want to duplicate this value - a weight index should be 0 or larger but it's optional to specify it

// general result for reading metrics - set this up to share for both metrics or comparisons (bounds are percents)
struct MetricAnalyzerInfo
{
  zscalar_t theValue;
  double upperValue;
  double lowerValue;
  bool bFoundUpperBound;
  bool bFoundLowerBound;
  std::string parameterDescription;
};

class MetricAnalyzer {
public:
  /// \brief Analyze metrics for a problem based on a range of tolerances
  ///
  /// @param metricsPlist parameter list defining tolerances
  /// @param problem the problem whose metrics are to be analyzed
  /// @param[out] msg_stream a std::ostringstream stream to return information from the analysis
  ///
  /// @return returns a boolean value indicated pass/failure.
  static bool analyzeMetrics(
    const RCP<const Zoltan2::EvaluatePartition <basic_id_t> > &metricObject,
    const ParameterList &metricsParameters, 
    std::ostringstream & msg_stream )
  {
    if (metricsParameters.numParams() == 0) {
      return true; // specification is that we do nothing - we may just be testing our status
    }

    bool bAllPassed = true;

    std::vector<MetricAnalyzerInfo> metricInfoSet;
    LoadMetricInfo(metricInfoSet, metricObject, metricsParameters);

    int countFailedMetricChecks = 0;
    for (auto metricInfo = metricInfoSet.begin();
         metricInfo != metricInfoSet.end(); ++metricInfo) {
      if (!MetricAnalyzer::executeMetricCheck(*metricInfo, msg_stream)) {
        ++countFailedMetricChecks;
      }
    }

    // this code prints a summary of all metric checks and indicates how many failed, if any did fail
    if(countFailedMetricChecks == 0) {
      msg_stream << metricsParameters.numParams() << " out of " << metricsParameters.numParams() << " metric checks"  << " PASSED." << std::endl;
    }
    else {
      msg_stream << countFailedMetricChecks << " out of " << metricsParameters.numParams() << " metric checks " << " FAILED." << std::endl;
      bAllPassed = false;
    }
    msg_stream << std::endl; // cosmetic spacer
    return bAllPassed;
  }

  static void LoadMetricInfo(
    std::vector<MetricAnalyzerInfo> & metricInfoSet,
    const RCP<const Zoltan2::EvaluatePartition <basic_id_t> > &metricObject,
    const ParameterList &metricsParameters) {

    // at this point we should be looking at a metricsPlist with the following 
    // format - note that weight is optional

    //      <ParameterList name="metriccheck1">
    //        <Parameter name="check" type="string" value="imbalance"/>
    //        <Parameter name="lower" type="double" value="0.99"/>
    //        <Parameter name="upper" type="double" value="1.4"/>
    //      </ParameterList>
    //      <ParameterList name="metriccheck2">
    //        <Parameter name="check" type="string" value="imbalance"/>
    //        <Parameter name="weight" type="int" value="0"/>
    //        <Parameter name="lower" type="double" value="0.99"/>
    //        <Parameter name="upper" type="double" value="1.4"/>
    //      </ParameterList>

    // first let's get a list of all the headings, so "metriccheck1", "metriccheck2" in this case
    // I've currently got this enforcing those names strictly to make sure formatting is correct
    // But really the headings could just be any unique names and are arbitrary
    int headingIndex = 1;

    for (auto iterateArbitraryHeadingNames = metricsParameters.begin(); 
         iterateArbitraryHeadingNames != metricsParameters.end();
         ++iterateArbitraryHeadingNames) {
      auto headingName = metricsParameters.name(iterateArbitraryHeadingNames);

      // we could be flexible on these headers but for now let's enforce it to get any convention inconsistencies cleaned up
      std::string expectedHeadingName = "metriccheck" + std::to_string(headingIndex);
      if( expectedHeadingName != headingName) {
        throw std::logic_error( "The parameter list expected to find a heading with name '" + expectedHeadingName + "' but instead found '" + headingName );
      }

      // get the parameters specific to the check we want to run
      const ParameterList & metricCheckParameters = metricsParameters.sublist(headingName);

      MetricAnalyzerInfo metricInfo = getMetricInfo(metricCheckParameters, metricObject);
      metricInfoSet.push_back(metricInfo);
      ++headingIndex;
    }
  }

private:

  static zscalar_t convertParameterChoicesToEvaluatePartitionAPICall( 
    const RCP<const Zoltan2::EvaluatePartition <basic_id_t> > &metricObject,
    std::string keyWord,
    int weightIndex,
    int selectedNormedSetting ) 
  {
    // this is going to need some consideration - how is the adapter scalar_t type to be properly handled?
    zscalar_t theValue = 0;

    if( weightIndex != UNDEFINED_PARAMETER_INT_INDEX && selectedNormedSetting != UNDEFINED_PARAMETER_INT_INDEX ) {
      throw std::logic_error( "Both parameters 'normed' and 'weight' were specified. They should never appear together." );
    }

    // this may not always be true - will have to find out how normed may be used late
    if( keyWord != API_STRING_getWeightImbalance && selectedNormedSetting != UNDEFINED_PARAMETER_INT_INDEX ) {
      throw std::logic_error( "'normed' was specified but this only has meaning for the 'imbalance' parameter." );
    }

    // Here I am enforcing a parallel usage to the way API calls exist in EvaluatePartition
    if (keyWord == API_STRING_getWeightImbalance) {
      if( weightIndex == UNDEFINED_PARAMETER_INT_INDEX ) {	// -1 is the code for optional (meaning it was not specified)
        if( selectedNormedSetting == 1 ) {
          theValue = metricObject->getNormedImbalance();
        }
        else {
          theValue = metricObject->getObjectCountImbalance(); // this will be index
        }
      }
      else {
        theValue = metricObject->getWeightImbalance(weightIndex); // this will get the proper index specified
      }
    }
    else if (keyWord == API_STRING_getTotalEdgeCuts) {
      if( weightIndex == UNDEFINED_PARAMETER_INT_INDEX ) {
        theValue = metricObject->getTotalEdgeCut();
      }
      else {
        theValue = metricObject->getTotalWeightEdgeCut(weightIndex);
      }
    }
    else if (keyWord == API_STRING_getMaxEdgeCuts) {
      if( weightIndex == UNDEFINED_PARAMETER_INT_INDEX ) {
        theValue = metricObject->getMaxEdgeCut();
      }
      else {
        theValue = metricObject->getMaxWeightEdgeCut(weightIndex);
      }
    }
    else {
      // we have found an invalid key word - throw an error
      throw std::logic_error( "The parameter '" + std::string(KEYWORD_PARAMETER_NAME) + "' was specified as '" + keyWord + "' which is not understood." );
    }

    return theValue;
  }

  static bool executeMetricCheck(
    const MetricAnalyzerInfo & metricInfo, 
    std::ostringstream &msg_stream) 
  {
    bool bDoesThisTestPass = true; // will set this false if a test fails
    if (metricInfo.bFoundUpperBound && metricInfo.bFoundLowerBound) {
      if (metricInfo.theValue < metricInfo.lowerValue || 
          metricInfo.theValue > metricInfo.upperValue) {
        msg_stream << "FAILED: " << metricInfo.parameterDescription 
                   << " value: " << metricInfo.theValue << " is not in range: "                    << metricInfo.lowerValue << " to " 
                   << metricInfo.upperValue << std::endl;
        bDoesThisTestPass = false;
      }
      else {
        msg_stream << "Success: " << metricInfo.parameterDescription  
                   << " value: " << metricInfo.theValue << " is in range: " 
                   << metricInfo.lowerValue << " to " 
                   << metricInfo.upperValue << std::endl;
      }
    }
    else if (metricInfo.bFoundUpperBound) {
      if (metricInfo.theValue > metricInfo.upperValue) {
        msg_stream << "FAILED: " << metricInfo.parameterDescription 
                   << " value: " << metricInfo.theValue << " is not below " 
                   << metricInfo.upperValue << std::endl;
        bDoesThisTestPass = false;
      }
      else {
        msg_stream << "Success: " << metricInfo.parameterDescription  
                   << " value: " << metricInfo.theValue << " is below: " 
                   << metricInfo.upperValue << std::endl;
      }
    }
    else if (metricInfo.bFoundLowerBound) {
      if (metricInfo.theValue < metricInfo.lowerValue) {
        msg_stream << "FAILED: " << metricInfo.parameterDescription 
                   << " value: " << metricInfo.theValue << " is not above " 
                   << metricInfo.lowerValue << std::endl;
        bDoesThisTestPass = false;
      }
      else {
        msg_stream << "Success: " << metricInfo.parameterDescription  
                   << " value: " << metricInfo.theValue << " is above: " 
                   << metricInfo.lowerValue << std::endl;
      }
    }
    return bDoesThisTestPass;
  }

  static MetricAnalyzerInfo getMetricInfo(  
    const ParameterList & metricCheckParameters,
    const RCP<const Zoltan2::EvaluatePartition <basic_id_t> > &metricObject)
  {
    MetricAnalyzerInfo result; // will fill these values
    for (auto iterateAllKeys = metricCheckParameters.begin();
         iterateAllKeys != metricCheckParameters.end(); ++iterateAllKeys) {
      auto checkName = metricCheckParameters.name(iterateAllKeys);
      if ( checkName != WEIGHT_PARAMETER_NAME &&
          checkName != KEYWORD_PARAMETER_NAME &&
          checkName != UPPER_PARAMETER_NAME &&
          checkName != LOWER_PARAMETER_NAME &&
          checkName != NORMED_PARAMETER_NAME ) {
        throw std::logic_error( "Key name: '" + checkName + "' is not understood." );
      }
    }

    // pick up the weight index - this parameter is optional so we check it first - this way we can communicate with EvaluatePartition properly in the next step
    int selectedWeightIndex = UNDEFINED_PARAMETER_INT_INDEX; // meaning not specified so default case
    if( metricCheckParameters.isParameter(WEIGHT_PARAMETER_NAME)) {
      selectedWeightIndex = metricCheckParameters.get<int>(WEIGHT_PARAMETER_NAME);
      if( selectedWeightIndex < 0 ) {
        throw std::logic_error( "Optional weight index was specified as: " + std::to_string(selectedWeightIndex) + "   Weight index must be 0 or positive." ); // I think that's the best I can do for error checking weight index right now - we may want to specify the cap when we know it
      }
    }

    // pick up the norm index - this parameter is optional so we check it first - this way we can communicate with EvaluatePartition properly in the next step
    int selectedNormedSetting = UNDEFINED_PARAMETER_INT_INDEX; // meaning not specified so default case
    if( metricCheckParameters.isParameter(NORMED_PARAMETER_NAME)) {
      bool bNormSetting = metricCheckParameters.get<bool>(NORMED_PARAMETER_NAME);
      selectedNormedSetting = bNormSetting ? 1 : 0;
      if( selectedNormedSetting != 0 && selectedNormedSetting != 1 ) {
        throw std::logic_error( "Optional normed parameter was specified as: " + std::to_string(selectedNormedSetting) + "   Normed parameter must be true or false." );
      }
    }

    // one of the parameters called "check" should define a string which is a keyword which correlates to an EvaluatePartition API all
    // this area needs some consideration - how and where shall we map and define the naming conventions
    if( metricCheckParameters.isParameter(KEYWORD_PARAMETER_NAME)) {
      std::string theKeyWord = metricCheckParameters.get<std::string>(KEYWORD_PARAMETER_NAME);

      // this is going to need some consideration - how is the adapter scalar_t type to be properly handled?
      result.theValue = convertParameterChoicesToEvaluatePartitionAPICall(metricObject, theKeyWord, selectedWeightIndex, selectedNormedSetting );

      // now we can obtain the upper and lower bounds for this test
      result.bFoundUpperBound = metricCheckParameters.isParameter(UPPER_PARAMETER_NAME);
      result.bFoundLowerBound = metricCheckParameters.isParameter(LOWER_PARAMETER_NAME);


      if (result.bFoundUpperBound) {
        result.upperValue = metricCheckParameters.get<double>(UPPER_PARAMETER_NAME);
      }
      if (result.bFoundLowerBound) {
        result.lowerValue = metricCheckParameters.get<double>(LOWER_PARAMETER_NAME);
      }

      result.parameterDescription = theKeyWord;
      if( selectedWeightIndex != UNDEFINED_PARAMETER_INT_INDEX ) {
        result.parameterDescription = result.parameterDescription + " (weight: " + std::to_string(selectedWeightIndex) + ")";
      }
      else if( selectedNormedSetting != UNDEFINED_PARAMETER_INT_INDEX ) {	// throw above would catch the case where both of these were set
        result.parameterDescription = result.parameterDescription + " (normed: " + ( ( selectedNormedSetting == 0 ) ? "false" : "true" ) + ")";
      }
    }
    return result;
  }
};


#endif //ZOLTAN2_METRIC_ANALYZER_HPP

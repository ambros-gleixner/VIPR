/*
*
*   Copyright (c) 2021 Fabian Frickenstein
*
*   Permission is hereby granted, free of charge, to any person obtaining a
*   copy of this software and associated documentation files (the "Software"),
*   to deal in the Software without restriction, including without limitation
*   the rights to use, copy, modify, merge, publish, distribute, sublicense,
*   and/or sell copies of the Software, and to permit persons to whom the
*   Software is furnished to do so, subject to the following conditions:
*
*   The above copyright notice and this permission notice shall be included in
*   all copies or substantial portions of the Software.
*
*   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
*   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
*   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
*   DEALINGS IN THE SOFTWARE.
*
*/

// Version control
#define VERSION_MAJOR 1
#define VERSION_MINOR 0

// Includes
#include <iostream>
#include <vector>
#include <map>
#include <limits>
#include "soplex.h"

using namespace std;
using namespace soplex;

// sparse boolean vectors
typedef map<int, bool> SVectorBool;

// Sparse vectors of rational numbers as maps
class SVectorRat : public map<int, Rational>
{
   public:
      void compactify() {  if( !_compact )
                           {
                              auto it = this->begin();
                              while( it != this->end() )
                              {
                                 if( it->second == 0 ) this->erase(it++);
                                 else ++it;
                              }
                              _compact = true;
                           }
                        }
      bool operator!=(SVectorRat &other);
      bool operator==(SVectorRat &other) { return !(*this != other);}

   private:
      bool _compact = false;
};

// Globals
ifstream certificateFile;
ofstream completedFile;
bool debugmode = false;
int numberOfVariables = 0; // number of variables
int numberOfIntegers = 0; // number of integers
int numberOfConstraints = 0; // number of constraints
int numberOfBoundedCons = 0; // number of bounded constraints
int numberOfSolutions = 0; // number of solutions
long numberOfDerivations = 0; // number of derivations
long currentDerivation = 0; // current derivation
DSVectorRational dummycol(0); // SoPlex placeholder column
vector<string> variables; // variable names
vector<bool> isInt; // integer variable indices
DSVectorRational ObjCoeff(0); // sparse vector of objective coefficients

vector<long> currentlyActiveDerivations; // tracks curr. active derivations (for completing lin)

vector<tuple<DSVectorRational, Rational, int>> constraints; // all constraints, including derived ones

struct lpIndex { long idx; bool isRowId; };
bool operator< ( lpIndex a, lpIndex b ) {
   return make_pair( a.idx, a.isRowId ) < make_pair( b.idx, b.isRowId );
}

typedef pair<long, long> certIndices;
map< lpIndex, certIndices> correspondingCertRow; // maps each row/var of the LP to its active con/asm
map< lpIndex, certIndices> originalCertRow; // maps each row/var of the LP to its original con/asm
VectorRational dualmultipliers(0);
VectorRational reducedCosts(0);
vector<tuple<Rational,Rational,long>> lowerBounds; // rational boundval, multiplier, long certindec
vector<tuple<Rational,Rational,long>> upperBounds; // rational boundval, multiplier, long certindex

map< long, lpIndex > correspondingLpData; // maps each certificate row to the corresponding lp Data

// Forward Declaration
void modifyFileName(string &path, const string &newExtension);
bool checkversion(string ver);
bool processVER();
bool processVAR(SoPlex &workinglp);
bool processINT();
bool processOBJ(SoPlex &workinglp);
bool processCON(SoPlex &workinglp);
bool processRTP();
bool processSOL();
bool processDER(SoPlex &workinglp);
bool getConstraints(SoPlex &workinglp, string &consense, Rational &rhs, int &activeConstraint);
bool derisasm();
bool derislin(SoPlex &workinglp, DSVectorRational &row, string &consense, Rational &rhs, string &label);
bool derisrnd();
bool derisuns();
bool derissol();
bool completeLin(SoPlex &workinglp, vector<long> &derToDelete, vector<long> &derToAdd, string &label);
bool printReasoningToCertificate(DVectorRational &dualmultipliers, DVectorRational &reducedcosts);

static void processGlobalBoundChange(Rational rhs, Rational boundmult, int varindex,
                                       long boundindex, int sense)
{
   char restline[256];
   string tmpstring;
   string lastword;

   if( certificateFile.peek() != '\n' )
   {
      certificateFile.getline(restline, 256);
      tmpstring = restline;
      istringstream tmpstream(tmpstring);

      while(tmpstream >> lastword)
      {
         if( lastword == "global" )
         {
            if( sense <= 0 ) // eq or ub
               upperBounds[varindex] = make_tuple(rhs, boundmult, boundindex);
            if( sense >= 0 ) // eq or lb
               lowerBounds[varindex] = make_tuple(rhs, boundmult, boundindex);
         }
      }
   }
}

// Main function
int main(int argc, char *argv[])
{
   int returnStatement = -1;

   if( argc < 2 || argc > 4)
   {
      cerr << "Usage: " << argv[0]
      << " <certificate filename> + SoPlex verbosity (optional) + Debugmode ON/OFF (optional)" << endl;
      return returnStatement;
   }

   certificateFile.open(argv[1]);

   if( certificateFile.fail() )
   {
      cerr << "Failed to open file " << argv[1] << endl;
      return returnStatement;
   }

   string path = argv[1];
   modifyFileName(path, "_complete.vipr");

   completedFile.open( path.c_str(), ios::out );

   if( completedFile.fail() )
   {
      cerr << "Failed to open file " << path << endl;
      return returnStatement;
   }

   // Set parameters for exact solving
   SoPlex baselp;
   baselp.setIntParam(SoPlex::READMODE, SoPlex::READMODE_RATIONAL);
   baselp.setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
   baselp.setIntParam(SoPlex::CHECKMODE, SoPlex::CHECKMODE_RATIONAL);
   baselp.setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
   baselp.setRealParam(SoPlex::FEASTOL, 0.0);
   baselp.setRealParam(SoPlex::OPTTOL, 0.0);

   if( argc == 2 )
   {
      baselp.setIntParam(SoPlex::VERBOSITY, 0);
      cout << "Verbosity level set to 0 (default)." << endl;
      cout << "Debugmode set to OFF (default)." << endl;
   }

   else if( argc == 3 )
   {

      if( isdigit(*argv[2]) )
      {
         int verbosity;
         verbosity = atoi(argv[2]);
         if( verbosity < 0 || verbosity > 5 )
         {
            baselp.setIntParam(SoPlex::VERBOSITY, 0);
            cerr << "Verbosity level outside range 0 to 5. Read " << verbosity << " instead." << endl;
            cout << "Continue with default settings (verbosity = 0)." << endl;
         }
         else
         {
            baselp.setIntParam(SoPlex::VERBOSITY, verbosity);
            cout << "Verbosity level set to " << verbosity << "." << endl;
         }
      }
      else
      {
         if( string(argv[2]) == "ON")
         {
            baselp.setIntParam(SoPlex::VERBOSITY, 0);
            debugmode = true;
            cout << "Verbosity level set to 0 (default)." << endl;
            cout << "Debugmode turned on." << endl;
         }
         else if( string(argv[2]) == "OFF")
         {
            baselp.setIntParam(SoPlex::VERBOSITY, 0);
            cout << "Verbosity level set to 0 (default)." << endl;
            cout << "Debugmode turned off." << endl;
         }
         else
         {
            baselp.setIntParam(SoPlex::VERBOSITY, 0);
            cout << "Unknown input for SoPlex verbosity (range 0 to 5 expected) "
            << "or debug settings (ON/OFF expected). Read "
            << string(argv[2]) << " instead" << endl;
            cout << "Verbosity level set to 0 (default)." << endl;
            cout << "Debugmode set to OFF (default)." << endl;
         }
      }
   }
   else if( argc == 4 )
   {
      // Set SoPlex verbosity
      int verbosity;
      if( !isdigit(*argv[2]) )
      {
         baselp.setIntParam(SoPlex::VERBOSITY, 0);
         cout << "Verbosity level is not a number. Read " << argv[2] << " instead." << endl;
         cout << "Continue with default settings (verbosity = 0)." << endl;
      }
      verbosity = atoi(argv[2]);
      if( verbosity < 0 || verbosity > 5 )
      {
         baselp.setIntParam(SoPlex::VERBOSITY, 0);
         cerr << "Verbosity level outside range 0 to 5. Read " << verbosity << " instead." << endl;
         cout << "Continue with default settings (verbosity = 0)." << endl;
      }
      else
      {
         baselp.setIntParam(SoPlex::VERBOSITY, verbosity);
         cout << "Verbosity level set to " << verbosity << "." << endl;
      }
      // Set Debugmode
      if( string(argv[3]) == "ON")
      {
         debugmode = true;
         cout << "Debugmode turned on." << endl;
      }
      else if( string(argv[3]) == "OFF")
      {
         cout << "Debugmode turned off." << endl;
      }
      else
      {
         cout << "Unknown input for debug settings (ON/OFF expected). Read "
         << string(argv[3]) << " instead" << endl;
         cout << "Continue with default setings (debugmode OFF)" << endl;
      }
   }



   double start_cpu_tm = clock();
   if( processVER() )
      if( processVAR(baselp) )
         if( processINT() )
            if( processOBJ(baselp) )
               if( processCON(baselp) )
                  if( processRTP() )
                     if( processSOL() )
                        if( processDER(baselp) )
                        {
                              cout << "Completion of File successful!" <<endl;
                              returnStatement = 0;
                              double cpu_dur = (clock() - start_cpu_tm)
                                               / (double)CLOCKS_PER_SEC;

                              cout << endl << "Completed in " << cpu_dur
                                   << " seconds (CPU)" << endl;
                        }
   return returnStatement;
}

// append "_complete.vipr" to filename
void modifyFileName(string &path, const string& newExtension)
{
   string::size_type position = path.find_last_of('.');
   if( position == string::npos )
      position = path.size();
   path.replace(position, newExtension.length(), newExtension);
}

// Version control for .vipr input file. Backward compatibility possible for minor versions
bool checkVersion(string version)
{
   bool returnStatement = false;

   size_t position = version.find(".");
   completedFile << " " + version;
   int major = atoi(version.substr(0, position).c_str());
   int minor = atoi(version.substr(position+1, version.length()-position).c_str());

   cout << "Certificate format version " << major << "." << minor << endl;

   if( (major == VERSION_MAJOR) && (minor <= VERSION_MINOR) )
   {
      returnStatement = true;
   }
   else
   {
      cerr << "Version unsupported" << endl;
   }

   return returnStatement;
}

// Check version and correct file format allows next process to start
// Error if the version is incompatible or not specified
bool processVER()
{
   bool returnStatement = false;
   string tmpStr;

   for( ;; )
   {
      certificateFile >> tmpStr;
      if( tmpStr == "VER" )
      {
         completedFile << tmpStr;
         certificateFile >> tmpStr;

         returnStatement = checkVersion(tmpStr);
break;
      }
      else if( tmpStr == "%" )
      {
         getline( certificateFile, tmpStr );
      }
      else
      {
         cerr << endl << "Comment or VER expected. Read instead "
                << tmpStr << endl;
break;
      }
   }

   return returnStatement;
}

// Processes number of variables
// Puts placeholder variables in SoPlex LP
// Error if nr of variables invalid or number of variables > specified variables or section missing
bool processVAR(SoPlex &workinglp)
{
   cout << endl << "Processing VAR section..." << endl;

   auto returnStatement = true;

   string section;

   certificateFile >> section;



   // Check section
   if( section != "VAR" )
   {
      cerr << "VAR expected.   Read instead " << section << endl;
   }
   else
   {
      completedFile << "\r\n" + section;
      certificateFile >> numberOfVariables; // number of variables

      if( certificateFile.fail() || numberOfVariables < 0 )
      {
         cerr << "Invalid number after VAR" << endl;
         returnStatement = false;
      }

      // Store variables
      else
      {
         completedFile << " " << numberOfVariables;
         isInt.resize( numberOfVariables );

         upperBounds.resize( numberOfVariables );
         lowerBounds.resize( numberOfVariables );
         fill(upperBounds.begin(), upperBounds.end(), make_tuple(infinity, 1, -1));
         fill(lowerBounds.begin(), lowerBounds.end(), make_tuple(-infinity, 1, -1));
         for( int i = 0; i < numberOfVariables; ++i )
         {
            string tmp;
            certificateFile >> tmp;
            if( certificateFile.fail() )
            {
               cerr << "Error reading variable for index " << i << endl;
               returnStatement = false;
         break;
            }
            completedFile <<"\r\n" + tmp;
            workinglp.addColRational( LPColRational( 1, dummycol, infinity, -infinity ) );
            variables.push_back( tmp );
            correspondingCertRow[{i, false}] = make_pair( -1, -1 );
            originalCertRow[{i, false}] = make_pair( -1, -1);
            isInt[i] = false;
         }
      }
   }

   return returnStatement;
}

// Processes number of integer variables and specifies their indices
// Produces vector of indices of integer variables
// Error if nr of integers invalid or nr of integers > specified integers or section missing
bool processINT()
{
   cout << endl << "Processing INT section..." << endl;

   bool returnStatement = false;

   string section;
   certificateFile >> section;

   if( section!= "INT" )
   {
      cerr << "INT expected. Read instead: " << section << endl;
      return returnStatement;
   }

   completedFile << "\r\n" + section;
   certificateFile >> numberOfIntegers;


   if( certificateFile.fail() )
   {
      cerr << "Failed to read number after INT" << endl;
      return returnStatement;
   }

   completedFile << " " <<numberOfIntegers << "\r\n";
   int index;

   for( int i = 0; i < numberOfIntegers; ++i )
   {
      certificateFile >> index;
      if( certificateFile.fail() )
      {
         cerr << "Error reading integer index " << i << endl;
         return returnStatement;
      }
      completedFile << index << " ";
      isInt[index] = true;
   }
   returnStatement = true;

   return returnStatement;
}

// Processes the sense of the objective function and coefficients for variables
// Stores sense and coefficients for possible future use
// Error if objective sense invalid (other than -1, 0, 1 for min, equality, max) or subroutine fails
bool processOBJ(SoPlex &workinglp)
{
   cout << endl << "Processing OBJ section..." << endl;

   bool returnStatement = false;
   string section, objectiveSense;
   int numberOfObjCoeff, idx;
   vector<Rational> values;
   vector<int> indices;
   VectorRational Objective(0); // full objective Vector
   Rational val;


   certificateFile >> section;

   if( section != "OBJ" )
   {
      cerr << "OBJ expected. Read instead: " << section << endl;
      return returnStatement;
   }
   completedFile << "\r\n" + section;
   certificateFile >> objectiveSense;
   completedFile << " " + objectiveSense;
   if( objectiveSense == "min" )
   {
      workinglp.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);
   }
   else if( objectiveSense == "max" )
   {
      workinglp.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MAXIMIZE);
   }
   else
   {
      cerr << "Invalid objective sense: " << objectiveSense << endl;
      return returnStatement;
   }

   Objective.reSize(numberOfVariables);
   Objective.reDim(numberOfVariables);

   certificateFile >> numberOfObjCoeff;
   completedFile << "\r\n" << numberOfObjCoeff << " ";

   values.reserve(numberOfObjCoeff);
   indices.reserve(numberOfObjCoeff);

   for( int j = 0; j < numberOfObjCoeff; ++j)
   {
      certificateFile >> idx >> val;
      completedFile << idx << " " << val << " ";
      values.push_back(val);
      indices.push_back(idx);
   }

   ObjCoeff.add(numberOfObjCoeff, indices.data(), values.data());
   Objective.assign(ObjCoeff);
   workinglp.changeObjRational(Objective);

   returnStatement = true;
   return returnStatement;
}

// Processes the constraint section of the file
// Stores constraints for future access
bool processCON(SoPlex &workinglp)
{
   cout << endl << "Processing CON section..." << endl;

   bool returnStatement = false;
   string section;
   string label, consense;
   Rational rhs;


   certificateFile >> section;
   if( section!= "CON" )
   {
      cerr << "CON expected. Read instead: " << section << endl;
      return returnStatement;
   }
   completedFile << "\n" + section;

   certificateFile >> numberOfConstraints >> numberOfBoundedCons;
   completedFile << " " << numberOfConstraints << " " << numberOfBoundedCons;

   for( int i = 0; i < numberOfConstraints; ++i)
   {
      certificateFile >> label >> consense >> rhs;
      completedFile << "\r\n" + label + "  " + consense + " " << rhs << "  ";
      currentDerivation++;
      returnStatement = getConstraints(workinglp, consense, rhs, i);
   }

   return returnStatement;
}

// Proecesses the relation to prove in order to complete file
bool processRTP()
{
   cout << endl << "Processing RTP section..." << endl;

   bool returnStatement = false;

   string section, relationToProveTypeStr, upperString, lowerString;

   certificateFile >> section;

   if( section != "RTP" )
   {
      cerr << "RTP expected.  Read instead " << section << endl;
   }
   else
   {
      completedFile << "\r\n" + section;
      certificateFile >> relationToProveTypeStr;
      if( relationToProveTypeStr == "infeas")
      {
         completedFile << " " + relationToProveTypeStr;
         returnStatement = true;
      }
      else if( relationToProveTypeStr != "range" )
      {
         cerr << "RTP: unrecognized verification type: " << relationToProveTypeStr << endl;
         return returnStatement;
      }
      else
      {
         completedFile << " " + relationToProveTypeStr;
         certificateFile >> lowerString >> upperString;
         completedFile << " " + lowerString + " " + upperString;
         returnStatement = true;
      }

   }

   return returnStatement;
}

// Processes the solutions to check in order to complete file
bool processSOL()
{
   cout << endl << "Processing SOL section... " << endl;
   bool returnStatement = false;
   string section, label;
   int numberOfVarSol = 0;
   int idx;
   Rational val;

   certificateFile >> section;

   if( section != "SOL" )
   {
      cerr << "SOL expected.   Read instead " << section << endl;
      return returnStatement;
   }

   completedFile << "\r\n" + section;
   certificateFile >> numberOfSolutions;

   if( certificateFile.fail() )
   {
      cerr << "Failed to read number after SOL" << endl;
      return returnStatement;
   }
   else if( numberOfSolutions < 0 )
   {
      cerr << "Invalid number after SOL: " << numberOfSolutions << endl;
      return returnStatement;
   }
   else
   {
      completedFile << " " << numberOfSolutions;

      for( int i = 0; i < numberOfSolutions; ++i )
      {
         certificateFile >> label >> numberOfVarSol;
         completedFile << "\r\n" << label << " " << numberOfVarSol;
         for( int i = 0; i < numberOfVarSol; ++i )
         {
            certificateFile >> idx >> val;

            completedFile << "  " << idx << " " << val;
         }
      }

      returnStatement = true;
   }
   return returnStatement;
}

static bool isEqual(DSVectorRational row1, DSVectorRational row2)
{
   if( row1.size() != row2.size() )
      return false;
   else
   {
      for( int i = 0; i < row1.size(); i++ )
      {
         if( row1.index(i) != row2.index(i))
            return false;
         if( row1[i] != row2[i] )
            return false;
      }
   }
   return true;
}

// write a row to the the completed File.
// This does not write the reasoning why this constraint is valid
bool printRowToCertificate(DSVectorRational row, string& sense, Rational rhs, string label)
{
   completedFile << label << " ";
   completedFile << sense;

   completedFile << " " << rhs << " ";

   if( isEqual(row, ObjCoeff) )
   {
      completedFile << "OBJ ";
   }
   else
   {
      completedFile << row.size();
      for( size_t i = 0; i < row.size(); i++ )
      {
         completedFile << " " << row.index(i) << " " << row[row.index(i)];
      }
   }

   return true;
}

// Processes Derivation section
// Complete derivations marked "incomplete" or "weak"
bool processDER(SoPlex &workinglp)
{
   cout << endl << "Processing DER section... " << endl;
   bool returnStatement = false;
   string section, numberOfCoefficients, label, consense, bracket, kind;
   int intOfCoefficients = 0;
   long idx, sense;
   Rational val, rhs;
   currentDerivation = numberOfConstraints;

   certificateFile >> section;

   if( section != "DER" )
   {
      cerr << "DER expected.   Read instead " << section << endl;
      return false;
   }

   completedFile << "\r\n" + section;

   certificateFile >> numberOfDerivations;
   completedFile << " " << numberOfDerivations;

   if( numberOfDerivations == 0 )
   {
      cout << "Number of derivations = 0. Nothing to complete." << endl;
      return true;
   }


   for( long i = 0; i < numberOfDerivations; ++i )
   {
      DSVectorRational row(0);
      vector<Rational> values;
      vector<int> indices;
      currentDerivation += 1;

      int n = numberOfDerivations + numberOfConstraints;
      bool global;

      certificateFile >> label >> consense >> rhs;
      completedFile << endl;

      if( debugmode == true )
      {
         cout << "completing constraint " << label << endl;
      }
      // get derived constraints
      certificateFile >> numberOfCoefficients;

      if( numberOfCoefficients == "OBJ" )
      {
         row = ObjCoeff;
         intOfCoefficients = row.size();
      }
      else
      {
         intOfCoefficients = atoi(numberOfCoefficients.c_str());
         values.reserve(intOfCoefficients);
         indices.reserve(intOfCoefficients);
         for( int j = 0; j < intOfCoefficients; ++j )
         {
            certificateFile >> idx >> val;
            values.push_back(val);
            indices.push_back(idx);
         }

         row.add(intOfCoefficients, indices.data(), values.data());
      }

      switch(consense[0])
      {
         case 'E':
            sense = 0; break;
         case 'L':
            sense = -1; break;
         case 'G':
            sense = 1; break;
         default:
            cerr << "wrong sense for constraints " << consense << endl;
            break;
      }
      constraints.push_back(make_tuple(row, rhs, sense));

      // obtain derivation kind
      certificateFile >> bracket >> kind;

      if( bracket != "{" )
      {
         cerr << "Expecting { but read instead " << bracket << endl;
         return false;
      }
      else
      {
         if( kind == "asm")
         {
            printRowToCertificate(row, consense, rhs, label);
            completedFile << " " + bracket + " " + kind;
            returnStatement = derisasm();
         }
         else if( kind == "lin")
         {
            printRowToCertificate(row, consense, rhs, label);
            completedFile << " " + bracket + " " + kind;
            returnStatement = derislin(workinglp, row, consense, rhs, label);
            if( !returnStatement )
               cerr << "Could not process constraint " << label << endl;
         }
         else if( kind == "rnd" )
         {
            printRowToCertificate(row, consense, rhs, label);
            completedFile << " " + bracket + " " + kind;
            returnStatement = derisrnd();
         }
         else if( kind == "uns" )
         {
            printRowToCertificate(row, consense, rhs, label);
            completedFile << " " + bracket + " " + kind;
            returnStatement = derisuns();
         }
         else if( kind == "sol")
         {
            printRowToCertificate(row, consense, rhs, label);
            completedFile << " " + bracket + " " + kind;
            returnStatement = derisasm();
         }
         else if( kind == "asm")
         {
            printRowToCertificate(row, consense, rhs, label);
            completedFile << " " + bracket + " " + kind;
            returnStatement = derissol();
         }
         else
         {
            cerr << "Unknown reason. Nothing to complete." << endl;
            returnStatement = false;
         }
         if( intOfCoefficients == 1 )
            processGlobalBoundChange(rhs, row[row.index(0)], idx, numberOfConstraints + i, sense);
         else
            certificateFile.ignore(numeric_limits<streamsize>::max(), '\n');
      }
      if( !returnStatement )
         break;
   }
   return returnStatement;
}

// Read Coefficients for any constraint
// Modify main LP
bool getConstraints(SoPlex &workinglp, string &consense, Rational &rhs, int &activeConstraint)
{
   bool returnStatement = true;
   string numberOfCoefficients,normalizedSense, actualsense;
   int intOfCoefficients = 0, sense;
   DSVectorRational row(0);
   vector<Rational> values;
   vector<int> indices;
   long idx, lastrow;
   Rational val, normalizedRhs;

   certificateFile >> numberOfCoefficients;
   completedFile << " " + numberOfCoefficients;


   if( numberOfCoefficients == "OBJ" )
   {
      row = ObjCoeff;
   }
   else
   {
      intOfCoefficients = atoi(numberOfCoefficients.c_str());
      values.reserve(intOfCoefficients);
      indices.reserve(intOfCoefficients);
      for( int j = 0; j < intOfCoefficients; ++j )
      {
         certificateFile >> idx >> val;
         completedFile << " " << idx << " " << val;
         values.push_back(val);
         indices.push_back(idx);
      }
      if( intOfCoefficients == 1)
      {

         // normalize bound constraints

         if( consense == "E")
            normalizedSense = "E";
         else if( consense == "L" )
            if( sign(val) <= 0 )
               normalizedSense = "G";
            else normalizedSense = "L";
         else
            if( sign(val) <= 0 )
               normalizedSense = "L";
            else normalizedSense = "G";

         normalizedRhs = rhs/val;

         // check if normalized bound is an improvement over existing bound

         if( normalizedSense == "E" )
         {
            if( get<0>(upperBounds[idx]) > normalizedRhs )
            {
               upperBounds[idx] = make_tuple(normalizedRhs, Rational(val), currentDerivation -1);
            }
            if( get<0>(lowerBounds[idx]) < normalizedRhs )
            {
               lowerBounds[idx] = make_tuple(normalizedRhs, Rational(val), currentDerivation -1);
            }
         }
         else if( normalizedSense == "L" && get<0>(upperBounds[idx]) > normalizedRhs )
         {
            upperBounds[idx] = make_tuple(normalizedRhs, Rational(val), currentDerivation -1);
         }
         else if( normalizedSense == "G" && get<0>(lowerBounds[idx]) < normalizedRhs )
         {
            lowerBounds[idx] = make_tuple(normalizedRhs, Rational(val), currentDerivation -1);
         }
      }

      row.add(intOfCoefficients, indices.data(), values.data());
   }

   if( consense == "E")
      {
         workinglp.addRowRational( LPRowRational( rhs, row, rhs) );
         sense = 0;
      }
   else if( consense == "L" )
      {
         workinglp.addRowRational( LPRowRational( -infinity, row, rhs) );
         sense = -1;
      }
   else if( consense == "G" )
      {
         workinglp.addRowRational( LPRowRational( rhs, row, infinity) );
         sense = 1;
      }
   else
      returnStatement = false;

   lastrow = workinglp.numRows();
   correspondingCertRow[{lastrow-1, true}] = make_pair(activeConstraint, activeConstraint);
   constraints.push_back(make_tuple(row, Rational(rhs), sense));

   return returnStatement;
}
// Case derivation is assumption, original line is taken over
bool derisasm()
{
   string bracket;
   long derHierarchy;

   certificateFile >> bracket;
   if( bracket != "}" )
   {
      cerr << "Expecting } but read instead" << bracket << endl;
      return  false;
   }
   else
   {
      completedFile << " " + bracket;
      certificateFile >> derHierarchy;
      completedFile << " " << derHierarchy;
      return true;
   }
}

// Reads multipliers for completing weak domination
static bool readMultipliers( int &sense, SVectorRat &mult )
{

   int k;
   bool returnStatement = true;

   mult.clear();

   certificateFile >> k;

   for( auto j = 0; j < k; ++j )
   {
      Rational a;
      int index;

      certificateFile >> index >> a;

      if( a == 0 ) continue; // ignore 0 multiplier

      mult[index] = a;

      if( sense == 0 )
      {
         sense = get<2>(constraints[index]) * a.sign();
      }
      else
      {
         int tmp = get<2>(constraints[index]) * a.sign();
         if( tmp != 0 && sense != tmp )
         {
            cerr << "Coefficient has wrong sign for index " << index << endl;
            returnStatement = false;
            goto TERMINATE;
         }
      }
   }

TERMINATE:
   return returnStatement;
}


// Reads linear combinations for completing weak domination
static bool readLinComb( int &sense, Rational &rhs, SVectorRat& coefficients, SVectorRat& mult,
                  int currentConstraintIndex, SVectorBool &assumptionList)
{
   bool returnStatement = true;

   if( !readMultipliers(sense, mult) )
   {
      returnStatement = false;
   }
   else
   {
      rhs = 0;
      coefficients.clear();
      assumptionList.clear();

      for( auto it = mult.begin(); it != mult.end(); ++it )
      {
         auto index = it->first;
         auto a = it->second;

         // auto myassumptionList = constraint[index].getassumptionList();

         // for( auto it2 = myassumptionList.begin(); it2 != myassumptionList.end(); ++it2 )
         //    assumptionList[it2->first] = true;

         auto &con = constraints[index];

         DSVectorRational c = get<0>(con);

         for( auto i = 0; i < c.size(); ++i )
         {
            (coefficients)[c.index(i)] += a * c[c.index(i)];
         }

         rhs += a * get<1>(con);
      }
   }

   return returnStatement;
}

// Complete "lin"-type derivations marked "weak"
static bool completeWeakDomination(DSVectorRational &row, string &consense, Rational &rhs)
{
   SVectorRat coefDer;
   SVectorRat multDer;
   Rational rhsDer;
   Rational correctedSide;
   SVectorBool asmlist;
   int senseDer = 0;
   int nbounds;
   bool success;
   tuple<Rational,Rational,long> tup;
   string bracket;
   map<int,tuple<long, Rational>> localLowerBoundsToUse;
   map<int,tuple<long, Rational>> localUpperBoundsToUse;
   Rational boundval;
   Rational boundfactor;
   long long boundindex;

   certificateFile >> bracket >> nbounds;
   for( size_t i = 0; i < nbounds; i++ )
   {
      int varIndex;
      long boundIndex;
      Rational val;
      string type;

      certificateFile >> type >> varIndex >> boundIndex >> val;
      if( type == "L" )
      {
         localLowerBoundsToUse[varIndex] = make_tuple(boundIndex, val);
      }
      else if( type == "U" )
      {
         localUpperBoundsToUse[varIndex] = make_tuple(boundIndex, val);
      }
      else
      {
         cerr << "type does not match L/U, but is instead " << type << endl;
         abort();
         return false;
      }
   }

   certificateFile >> bracket;

   switch(consense[0])
   {
      case 'E':
         senseDer = 0; break;
      case 'L':
         senseDer = -1; break;
      case 'G':
         senseDer = 1; break;
      default:
         cerr << "wrong sense for constraints " << consense << endl;
         break;
   }

   if( !readLinComb(senseDer, rhsDer, coefDer, multDer, 0, asmlist) )
   {
      certificateFile.ignore(numeric_limits<streamsize>::max(), '\n');
      return false;
   }

   correctedSide = rhsDer;

   for( auto it = coefDer.begin() ; it != coefDer.end() ; ++it )
   {
      int idx = it->first;
      Rational derivedVal = it->second;
      Rational valToDerive = row[idx];

      coefDer[idx] = valToDerive;

      if( derivedVal == valToDerive )
         continue;
      else
      {
         // multiplier for the bound
         Rational boundmult = valToDerive - derivedVal;
         bool islower;

         // con is <= -> need to use ub for positive, lb for negative boundmult
         if(consense == "L")
            islower = (boundmult <= 0);
         // con is >= -> need to use ub for negative, lb for positive boundmult
         else if(consense == "G")
            islower = (boundmult >= 0);
         // cons is == -> this can currently not be handled, would need to be split in two parts
         else if(consense == "E")
         {
            cerr << "  cannot complete weak dominated equality constraints" << endl;
            certificateFile.ignore(numeric_limits<streamsize>::max(), '\n');
            return false;
         }

         // get the correct bound index and value
         if( islower && (localLowerBoundsToUse.find(idx) != localLowerBoundsToUse.end()) )
         {
            boundindex = get<0>(localLowerBoundsToUse[idx]);
            boundval = get<1>(localLowerBoundsToUse[idx]);
            boundfactor = 1;
         }
         else if( (!islower) && (localUpperBoundsToUse.find(idx) != localUpperBoundsToUse.end()) )
         {
            boundindex = get<0>(localUpperBoundsToUse[idx]);
            boundval = get<1>(localUpperBoundsToUse[idx]);
            boundfactor = 1;
         }
         else
         {
            tup = islower ? lowerBounds[idx] : upperBounds[idx];
            boundfactor = get<1>(tup);
            boundval = get<0>(tup);
            boundindex = get<2>(tup);
         }

         multDer[boundindex] += boundmult / boundfactor;
         correctedSide += boundmult * boundval;

         if( debugmode == true)
         {
            if( boundmult != 0 )
            {
               cout << "    correcting variable " << variables[idx] << " (idx " << idx << ") by " <<
                        boundmult << " (" << static_cast<double>(boundmult) << ") using";
               if( islower )
               {
                  if( localLowerBoundsToUse.count(idx) )
                     cout << " local ";
                  cout << " lower bound ";
               }
               else
               {
                  if( localUpperBoundsToUse.count(idx) )
                     cout << " local ";
                  cout << " upper bound ";
               }
               cout << boundval << " (" << static_cast<double>(boundval) << ")" << endl;
            }
         }
      }
   }
   // now go the other way
   for( int i = 0; i < row.size(); ++i )
   {
      int idx = row.index(i);
      Rational derivedVal = coefDer[idx];
      Rational valToDerive = row[idx];

      if( derivedVal == valToDerive )
         continue;
      else
      {
         // multiplier for the bound
         Rational boundmult = valToDerive - derivedVal;
         bool islower;

         // con is <= -> need to use ub for positive, lb for negative boundmult
         if(consense == "L")
            islower = (boundmult <= 0);
         // con is >= -> need to use ub for negative, lb for positive boundmult
         else if(consense == "G")
            islower = (boundmult >= 0);
         // cons is == -> this can currently not be handled, would need to be split in two parts
         else if(consense == "E")
         {
            cerr << "  cannot complete weak dominated equality constraints" << endl;
            certificateFile.ignore(numeric_limits<streamsize>::max(), '\n');
            return false;
         }

         // get the correct bound index and value
         if( islower && (localLowerBoundsToUse.find(idx) != localLowerBoundsToUse.end()) )
         {
            boundindex = get<0>(localLowerBoundsToUse[idx]);
            boundval = get<1>(localLowerBoundsToUse[idx]);
            boundfactor = 1;
         }
         else if( (!islower) && (localUpperBoundsToUse.find(idx) != localUpperBoundsToUse.end()) )
         {
            boundindex = get<0>(localUpperBoundsToUse[idx]);
            boundval = get<1>(localUpperBoundsToUse[idx]);
            boundfactor = 1;
         }
         else
         {
            tup = islower ? lowerBounds[idx] : upperBounds[idx];
            boundfactor = get<1>(tup);
            boundval = get<0>(tup);
            boundindex = get<2>(tup);
         }

         multDer[boundindex] += boundmult / boundfactor;
         correctedSide += boundmult * boundval;

         if( debugmode == true )
         {
            if( boundmult != 0 )
            {
               cout << "    correcting variable " << variables[idx] << " (idx " << idx << ") by "
                    << boundmult << " (" << static_cast<double>(boundmult) << ") using";
               if( islower )
               {
                  if( localLowerBoundsToUse.count(idx) )
                     cout << " local ";
                  cout << " lower bound ";
               }
               else
               {
                  if( localUpperBoundsToUse.count(idx) )
                     cout << " local ";
                  cout << " upper bound ";
               }
               cout << boundval << " (" << static_cast<double>(boundval) << ")" << endl;
            }
         }
      }
   }

   if( debugmode == true )
   {
      cout.precision(numeric_limits<double>::max_digits10);
      cout << "  exact rhs before correction: " << rhsDer << " ("
           <<    static_cast<double>(rhsDer) << ")" << endl;

      cout << "  exact rhs after  correction: " << correctedSide << " ("
           <<    static_cast<double>(correctedSide) << ")" << endl;

      cout << "  rhs that was printed       : " << rhs << " ("
           <<    static_cast<double>(rhs) << ")" << endl;
   }

   // first case: < and the side is larger, second case: > and side is smaller
   if( (senseDer == -1 && correctedSide > rhs) || (senseDer == 1 && correctedSide < rhs) )
   {
      cerr.precision(numeric_limits<double>::max_digits10);
      cerr << "Constraint does not dominate original one." << endl << "  Corrected Side is "
           << correctedSide << "(" << static_cast<double>(correctedSide) << ")" << endl
           <<  "  Original rhs is " << rhs << "(" << static_cast<double>(rhs) << ")" << endl;

      cerr << "  difference: " << static_cast<double>(correctedSide - rhs) << endl;
      success = false;
   }
   else
      success = true;

   completedFile << " " << multDer.size();
   for(auto it = multDer.begin(); it != multDer.end(); it++)
   {
      completedFile << " " << it->first << " " << it->second;
   }

   completedFile << " } -1";

   return success;
}

// Case derivation is "lin"
// Completes cases "incomplete" or "weak"
bool derislin(SoPlex &workinglp, DSVectorRational &row, string &consense, Rational &rhs, string &label)
{

   string numberOfCoefficients, tmp, bracket;
   long intOfCoefficients, idx, derhir;

   vector<long> newActiveDerivations;
   vector<long> toDeleteDerivations;
   vector<long> toAddDerivations;

   Rational val;
   DSVectorRational reasoningRow(0);

   certificateFile >> numberOfCoefficients;

   // Prepares completion of "incomplete" reasoning by locally modifying the LP
   if( numberOfCoefficients == "incomplete" )
   {
      VectorRational newObjective(0);
      newObjective.reSize(numberOfVariables);
      newObjective.reDim(numberOfVariables);
      newObjective = row;
      workinglp.changeObjRational(newObjective);

      if( consense == "G" or consense == "E" )
         workinglp.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);
      else if( consense == "L")
         workinglp.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MAXIMIZE);
      else
         cerr << "Invalid sense: " << consense << endl;

      certificateFile >> tmp;

      while( tmp != "}" )
      {
         newActiveDerivations.push_back(stol(tmp));

         certificateFile >> tmp;
      }

      sort(currentlyActiveDerivations.begin(), currentlyActiveDerivations.end(), less<int>());
      sort(newActiveDerivations.begin(), newActiveDerivations.end(), less<int>());

      set_difference( currentlyActiveDerivations.begin(), currentlyActiveDerivations.end(),
                        newActiveDerivations.begin(), newActiveDerivations.end(),
                        inserter( toDeleteDerivations, toDeleteDerivations.begin() ) );

      set_difference( newActiveDerivations.begin(), newActiveDerivations.end(),
                        currentlyActiveDerivations.begin(), currentlyActiveDerivations.end(),
                        inserter( toAddDerivations, toAddDerivations.begin() ) );

      toDeleteDerivations = currentlyActiveDerivations;
      toAddDerivations = newActiveDerivations;
      currentlyActiveDerivations = newActiveDerivations;
      return completeLin(workinglp, toDeleteDerivations, toAddDerivations, label);
   }

   else if( numberOfCoefficients == "weak")
      return completeWeakDomination(row, consense, rhs);

   else
   {
      intOfCoefficients = atoi(numberOfCoefficients.c_str());
      if( !isdigit(numberOfCoefficients[0]) )
      {
         cerr << "number of coefficients is not a number " << numberOfCoefficients << endl;
         certificateFile.ignore(numeric_limits<streamsize>::max(), '\n');
         return false;
      }
      completedFile << " " << intOfCoefficients;

      for( int i = 0; i < intOfCoefficients; ++i)
      {
         certificateFile >> idx >> val;
         completedFile << " " << idx << " " << val;
      }

      certificateFile >> bracket;
      if( bracket != "}" )
      {
         cerr << "Expecting } but read instead" << bracket << endl;
         return  false;
      }

      completedFile << " " << bracket;
      certificateFile >> derhir;
      completedFile << " " << derhir;
      return true;

   }
}


bool derisrnd()
{
   long numberOfCoefficients, idx, dernr;
   Rational val;
   string bracket;


   certificateFile >> numberOfCoefficients;
   completedFile << " " << numberOfCoefficients;

   for( int i = 0; i < numberOfCoefficients; ++i )
   {
      certificateFile >> idx >> val;
      completedFile << " " << idx << " " << val;
   }

   certificateFile >> bracket;
   if( bracket != "}" )
   {
      cerr << "Expecting } but read instead" << bracket << endl;
      return  false;
   }
   else
   {
      completedFile << " " + bracket;
      certificateFile >> dernr;
      completedFile << " " << dernr;
      return true;
   }
}


bool derisuns()
{
   long i1, l1, i2, l2, derHierarchy;
   string bracket;

   certificateFile >> i1 >> l1 >> i2 >> l2;
   completedFile << " " << i1 << " " << l1 <<
                    " " << i2 << " " << l2;

   certificateFile >> bracket;
   if( bracket != "}" )
   {
      cerr << "Expecting } but read instead" << bracket << endl;
      return  false;
   }
   else
   {
      completedFile << " " + bracket;
      certificateFile >> derHierarchy;
      completedFile << " " << derHierarchy;
      return true;
   }
}

bool derissol()
{

   string bracket;
   long derHierarchy;

   certificateFile >> bracket;
   if( bracket != "}" )
   {
      cerr << "Expecting } but read instead" << bracket << endl;
      return  false;
   }
   else
   {
      completedFile << " " + bracket;
      certificateFile >> derHierarchy;
      completedFile << " " << derHierarchy;
      return true;
   }
}

// Completes "incomplete" derivations
// Locally modifies LP and solves
bool completeLin(SoPlex &workinglp, vector<long> &derToDelete, vector<long> &derToAdd, string &label)
{
   string tmp;
   long numrows, derHierarchy;
   DSVectorRational row(0);
   int consense, normalizedSense;
   Rational rhs;
   lpIndex lpData;
   Rational normalizedRhs;
   vector<long> rowsToDelete;

   DVectorRational dualmultipliers(0);
   DVectorRational reducedcosts(0);

   vector<long>::iterator derIterator;

   // Delete rows and reset bounds from LP which are not used in the current completion attempt
   for( derIterator = derToDelete.begin(); derIterator != derToDelete.end(); derIterator++)
   {
      lpData = correspondingLpData[*derIterator];

      if( lpData.isRowId == true)
      {

         rowsToDelete.push_back(lpData.idx);
         correspondingCertRow.erase({lpData.idx, true});
         correspondingLpData.erase(*derIterator);
      }
      else
      {
         Rational originalUpperBound = get<0>(upperBounds[lpData.idx]);
         Rational originalLowerBound = get<0>(lowerBounds[lpData.idx]);

         workinglp.changeUpperRational( lpData.idx, infinity );
         workinglp.changeLowerRational( lpData.idx, -infinity );

         correspondingCertRow[{lpData.idx, false}] = originalCertRow[{lpData.idx, false}];
         correspondingLpData.erase(*derIterator);
      }
   }

   sort(rowsToDelete.begin(), rowsToDelete.end(), greater<int>());
   for( auto rowToDelete = rowsToDelete.begin(); rowToDelete != rowsToDelete.end(); rowToDelete++ )
   {
      workinglp.removeRowRational(*rowToDelete);
   }

   // Update LP with new derivations to be used
   for( derIterator = derToAdd.begin(); derIterator != derToAdd.end(); derIterator++ )
   {
      auto &missingCon = constraints[*derIterator];
      auto ncurrent = workinglp.numRowsRational();
      row = get<0>(missingCon);
      consense = get<2>(missingCon);
      rhs = get<1>(missingCon);

      if( row.size() == 1 && false )
      {
         // normalize bound constraints

         if( consense == 0)
            normalizedSense = 0;
         else if( consense == -1 )
            if( sign(row.value(0)) <= 0 )
               normalizedSense = 1;
            else normalizedSense = -1;
         else
            if( sign(row.value(0)) <= 0 )
               normalizedSense = -1;
            else normalizedSense = 1;

         normalizedRhs = rhs/row.value(0);

         if( normalizedSense == 0 )
         {
            if( workinglp.upperRational(row.index(0)) >= normalizedRhs )
            {
               if( workinglp.lowerRational(row.index(0)) > normalizedRhs )
               {
                  workinglp.addRowRational( LPRowRational( -infinity, row, rhs ) );
                  correspondingCertRow[{ workinglp.numRows()-1, true }] = make_pair( *derIterator, *derIterator );
                  correspondingLpData[*derIterator] = {workinglp.numRows()-1 , true};
                  assert(ncurrent + 1 == workinglp.numRowsRational());
               }

               else
               {
                  workinglp.changeUpperRational( row.index(0), normalizedRhs );
                  correspondingCertRow[{ row.index(0), false }] = make_pair( *derIterator, *derIterator );
                  correspondingLpData[*derIterator] = {row.index(0), false};
               }
            }
            else
            {
               cerr << "Error: in derivation " << label << ". New bound is no improvement.\n";
               return false;
            }

            if( workinglp.lowerRational(row.index(0)) <= normalizedRhs )
            {
               if( workinglp.upperRational(row.index(0)) < normalizedRhs )
               {
                  workinglp.addRowRational( LPRowRational( rhs, row, infinity ) );
                  correspondingCertRow[{ workinglp.numRows()-1, true }] = make_pair( *derIterator, *derIterator );
                  correspondingLpData[*derIterator] = {workinglp.numRows()-1 , true};
                  assert(ncurrent + 1 == workinglp.numRowsRational());
               }
               else
               {
                  workinglp.changeLowerRational( row.index(0), normalizedRhs );
                  correspondingCertRow[{ row.index(0), false }] = make_pair( *derIterator, *derIterator );
                  correspondingLpData[*derIterator] = {row.index(0), false};
               }
            }
            else
            {
               cerr << "Error: in derivation " << label << ". New bound is no improvement.\n";
               return false;
            }
         }



         else if( normalizedSense == -1 && workinglp.upperRational(row.index(0)) >= normalizedRhs )
         {
            if( workinglp.lowerRational(row.index(0)) > normalizedRhs )
            {
               workinglp.addRowRational( LPRowRational( -infinity, row, rhs ) );
               correspondingCertRow[{ workinglp.numRows()-1, true }] = make_pair( *derIterator, *derIterator );
               correspondingLpData[*derIterator] = {workinglp.numRows()-1 , true};
               assert(ncurrent + 1 == workinglp.numRowsRational());
            }

            else
            {
               workinglp.changeUpperRational( row.index(0), normalizedRhs );
               correspondingCertRow[{ row.index(0), false }] = make_pair( *derIterator, *derIterator );
               correspondingLpData[*derIterator] = {row.index(0), false};
            }
         }


         else if( normalizedSense == 1 && workinglp.lowerRational(row.index(0)) <= normalizedRhs )
         {
            if( workinglp.upperRational(row.index(0)) < normalizedRhs )
            {
               workinglp.addRowRational( LPRowRational( rhs, row, infinity ) );
               correspondingCertRow[{ workinglp.numRows()-1, true }] = make_pair( *derIterator, *derIterator );
               correspondingLpData[*derIterator] = {workinglp.numRows()-1 , true};
               assert(ncurrent + 1 == workinglp.numRowsRational());
            }
            else
            {
               workinglp.changeLowerRational( row.index(0), normalizedRhs );
               correspondingCertRow[{ row.index(0), false }] = make_pair( *derIterator, *derIterator );
               correspondingLpData[*derIterator] = {row.index(0), false};
            }
         }
         else
         {
            cerr << "Error: in derivation " << label << ". New bound is no improvement.";
            return false;
         }
      }
      else {
         if( consense == 0 )
               workinglp.addRowRational( LPRowRational( rhs, row, rhs ) );
         else if( consense == -1 )
               workinglp.addRowRational( LPRowRational( -infinity, row, rhs ) );
         else if( consense == 1 )
               workinglp.addRowRational( LPRowRational( rhs, row, infinity ) );
         else
            return false;

         correspondingCertRow[{ workinglp.numRows()-1, true }] = make_pair( *derIterator, *derIterator );
         correspondingLpData[*derIterator] = {workinglp.numRows()-1 , true};
         assert(ncurrent + 1 == workinglp.numRowsRational());
      }
   }



   SPxSolver::Status stat;
   numrows = workinglp.numRows();
   dualmultipliers.reDim(numrows);
   reducedcosts.reDim(numberOfVariables);

   stat = workinglp.optimize();

   if( stat == SPxSolver::OPTIMAL )
   {
      workinglp.getDualRational(dualmultipliers);
      workinglp.getRedCostRational(reducedcosts);

      printReasoningToCertificate(dualmultipliers, reducedcosts);
      completedFile << " }";
      certificateFile >> derHierarchy;
      completedFile << " " << derHierarchy;

   }
   else if( stat == SPxSolver::INFEASIBLE )
   {
      workinglp.getDualFarkasRational(dualmultipliers);
      workinglp.getRedCostRational(reducedcosts);

      printReasoningToCertificate(dualmultipliers, reducedcosts);
      completedFile << "}";
      certificateFile >> derHierarchy;
      completedFile << " " << derHierarchy;

   }
   else
   {
      cerr << "Warning: Completion attempt of Derivation "<< label <<" returned with status " << stat << ".\n";
      cerr << "Skip and continue completion of certificate.\n";
      completedFile << " incomplete";
      for( auto i: currentlyActiveDerivations )
         completedFile << " " << i;

   }
   return true;
}

bool printReasoningToCertificate(DVectorRational &dualmultipliers, DVectorRational &reducedcosts)
{
   DSVectorRational reasoningRow(0);
   long certIndex;
   Rational correctionFactor;

   for( int i= 0; i < reducedcosts.dim(); ++i )
   {
      if( sign(reducedcosts[i]) < 0 )
         certIndex = correspondingCertRow[{i, false}].second;
      else
         certIndex = correspondingCertRow[{i, false}].first;

      reasoningRow.add(certIndex, reducedcosts[i] );
   }

   for( int i = 0; i < dualmultipliers.dim(); ++i )
   {

      certIndex = correspondingCertRow[{i, true}].first;
      if( get<0>(constraints[certIndex]).dim() == 1 )
         correctionFactor = get<0>(constraints[certIndex]).value(0);
      else
         correctionFactor = 1;

      reasoningRow.add(certIndex, dualmultipliers[i] * correctionFactor);

   }
   completedFile << " " << reasoningRow.size();
   reasoningRow.sort();

   for( size_t i = 0; i < reasoningRow.size(); i++ )
      {
         completedFile << " " << reasoningRow.index(i) << " " << reasoningRow[reasoningRow.index(i)];
      }

   return true;
}
/*
*
*   Copyright (c) 2016 Kevin K. H. Cheung
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

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <assert.h>
#include <gmpxx.h>
#include <cstdio>
#include <ctime>
#include <chrono>

#define VERSION_MAJOR 1
#define VERSION_MINOR 0

using namespace std;

// Sparse vectors of rational numbers as maps
class SVectorGMP : public map<int, mpq_class>
{
   public:
      void compactify() {  if (!_compact) 
                           {
                              auto it = this->begin(); 
                              while (it != this->end())
                              {
                                 if (it->second == 0) this->erase(it++); 
                                 else ++it;
                              }
                              _compact = true;
                           }
                        }
      bool operator!=( SVectorGMP &other );
      bool operator==( SVectorGMP &other ) { return !(*this != other);}
   private:
      bool _compact = false;
};

typedef map<int, bool> SVectorBool;


// The type of derivation used to derive a constraint
enum DerType
{
   ASM,   // assumption
   LIN,   // simple implication
   RND,    // simple implication with integer rounding; i.e. a CG cut
   UNS,   // unsplit operation
   UNKNOWN 
};

enum RtpType
{
   INFEAS,
   RANGE
};


class Constraint;
class Constraint
{
   public:
      Constraint() {}

      Constraint( const string label, const int sense, const mpq_class rhs, 
                  const SVectorGMP coef, 
                  const bool isAsmCon, const SVectorBool asmList) : 
                  _label( label), _sense( sense ), _rhs (rhs), 
                  _coef(coef), _isAsm(isAsmCon), 
                  _asmList( asmList ) { _coef.compactify(); _trashed = false;
                                        _falsehood = _isFalsehood(); }


      bool round();

      mpq_class getRhs() { return _rhs; }
      mpq_class getCoef( const int idx ) { if (_coef.find(idx)!=_coef.end()) 
                                              return _coef[idx];
                                           else return mpq_class(0); }
      SVectorGMP coefSVec() { return _coef; }

      int getSense() { return _sense; }

      bool isAsm() { return _isAsm; }

      bool isFalsehood() { return _falsehood; } 
                  // true iff the constraint is a contradiction like 0 >= 1

      bool isTautology();
                  // true iff the constraint is a contradiction like 0 >= 1

      bool hasAsm( const int idx ) { return (_asmList.find(idx) != _asmList.end()); }

      void setAsmList( const SVectorBool asmList ) { _asmList = asmList; }
      SVectorBool getAsmList() { return _asmList; }

      bool dominates( Constraint &other );
      void print();

      void trash() { _trashed = true; _coef.clear(); }
      bool isTrashed() { return _trashed; }

      string label() { return _label; }

      void setMaxRefIdx( int refIdx ) { _refIdx = refIdx; }
      int getMaxRefIdx( ) { return _refIdx; }

   private:
      string _label;
      int _sense;
      mpq_class _rhs;
      SVectorGMP _coef;
      int _refIdx = -1;
      bool _isAsm;
      SVectorBool _asmList; // constraint index list that are assumptions
      bool _falsehood;

      bool _isFalsehood();
      bool _trashed;
};


// Globals
const SVectorBool emptyList;

int numVar = 0;
int numCon = 0;
int numBnd = 0;
int numDer = 0;
int numSol = 0;
vector<bool> isInt; // integer variable indices
vector<string> var; // variable names
vector<Constraint> constraint; // all the constraints, including derived ones
vector<SVectorGMP> solution; // all the solutions for checking feasibility
ifstream pf;   // certificate file stream

RtpType rtpType;
mpq_class bestObjVal; // best objective function value of specified solutions
mpq_class lowerBound; // lower bound for optimal value to be checked
mpq_class upperBound; // upper bound for optimal value to be checked
string lowerStr, upperStr;
bool isMin; // is minimization problem
bool checkLower; // true iff need to verify lower bound
bool checkUpper; // true iff need to verify upper bound
Constraint rtp; // constraint to be derived in the case of bound checking
SVectorGMP objCoef;

bool checkVersion(string ver);
bool processVER();
bool processVAR();
bool processINT();
bool processOBJ();
bool processCON();
bool processRTP();
bool processSOL();
bool processDER();

bool readMultipliers( int &sense, SVectorGMP &mult );
bool readConstraintCoefs( SVectorGMP &v );
bool readConstraint( string &label, int &sense, mpq_class &rhs, 
                               SVectorGMP &coef );
inline mpq_class floor( const mpq_class &q );
inline mpq_class ceil( const mpq_class &q );
bool isInteger( const mpq_class &q );

mpq_class scalarProd( const SVectorGMP &u, const SVectorGMP &v );
bool canUnsplit( Constraint &toDer, const int con1, const int a1, 
                 const int con2, const int a2, SVectorBool &asmList );

bool readLinComb( int &sense, mpq_class &rhs, SVectorGMP &coef, int currConIdx, 
                  SVectorBool &amsList ); 

int main( int argc, char *argv[] )
{

   int rs = -1;

   if (argc != 2)
   {
      cerr << "Usage: " << argv[0] << " <certificate filename>\n";
      return rs;
   }

   pf.open( argv[1] );

   if (pf.fail())
   {
      cerr << "Failed to open file " << argv[1] << endl;
      return rs;
   }

   double start_cpu_tm = clock();
   if( processVER() )
      if( processVAR() )
         if( processINT() )
            if( processOBJ() )
               if( processCON() )
                  if( processRTP() )
                     if( processSOL())
                        if( processDER() ) {
                           rs = 0;
                           double cpu_dur = (clock() - start_cpu_tm ) 
                                            / (double)CLOCKS_PER_SEC;

                           cout << endl << "Completed in " << cpu_dur 
                                << " seconds (CPU time)." << endl;
                        }

   return rs;
}

bool processRTP()
{

   cout << endl << "Processing RTP section..." << endl;

   bool rs = false;

   string section;

   pf >> section;
   
   if( section != "RTP" )
   {
      cerr << "RTP expected.   Read instead " << section << endl;
   }
   else
   {
      string rtpTypeStr;

      pf >> rtpTypeStr;

      if( rtpTypeStr == "infeas" )
      {
         rtpType = RtpType::INFEAS;

         cout << "Need to verify infeasibility. " << endl;
      }
      else if ( rtpTypeStr != "range" )
      {
         cerr << "RTP: unrecognized verification type: " << rtpTypeStr << endl;
         goto TERMINATE;
      }
      else
      {
         rtpType = RtpType::RANGE;

         checkLower = checkUpper = false;

         pf >> lowerStr >> upperStr;

         if( lowerStr != "-inf" )
         {
            checkLower = true;
            lowerBound = mpq_class( lowerStr );
         }

         if( upperStr != "inf" )
         {
            checkUpper = true;
            upperBound = mpq_class( upperStr );
         }

         if( checkLower && checkUpper && (lowerBound > upperBound) )
         {
            cerr << "RTP: invalid bounds. " << endl;
            goto TERMINATE;
         }       

         if ( isMin && checkLower )
         {
            rtp = Constraint("rtp", 1, lowerBound, objCoef, false, emptyList);
         }
         else if ( !isMin && checkUpper )
         {
            rtp = Constraint("rtp", -1, upperBound, objCoef, false, emptyList);
         }
         else
         {
            goto TERMINATE;
         }
   
         cout << "Need to verify optimal value range "
                << (lowerStr == "-inf" ? "(" : "[")
                << lowerStr << ", " << upperStr
                << (upperStr == "inf" ? ")" : "]")
                << "." << endl;
   
      }

      rs = true;
   }

TERMINATE:
   return rs;

}


bool checkVersion( string ver )
{
   bool rs = false;

   size_t pos = ver.find( "." );

   int major = atoi( ver.substr( 0, pos ).c_str() );
   int minor = atoi( ver.substr( pos+1, ver.length()-pos ).c_str() );

   cout << "Certificate format version " << major << "." << minor << endl;

   if ( (major == VERSION_MAJOR) && (minor <= VERSION_MINOR ) )
   {
      rs = true;
   }
   else
   {
      cerr << "Version unsupported" << endl;
   }

   return rs;
}


bool processVER()
{
   bool rs = false;
   string tmpStr;

   for(;;)
   {
      pf >> tmpStr;
      if( tmpStr == "VER" )
      {
         pf >> tmpStr;
         rs = checkVersion( tmpStr );
break;
      }
      else if( tmpStr == "%" )
      {
         getline( pf, tmpStr );
      }
      else
      {
         cerr << endl << "Comment or VER expected. Read instead " 
                << tmpStr << endl;
break;
      }
   }

   return rs;
}


bool processVAR()
{
   
   cout << endl << "Processing VAR section..." << endl;

   auto rs = true;

   string section;

   pf >> section;

   if ( section != "VAR" )
   {
      cerr << "VAR expected.   Read instead " << section << endl;
   }
   else
   {
      pf >> numVar; // number of variables
   
      if( pf.fail() || numVar < 0 )
      {
         cerr << "Invalid number after VAR" << endl;
         rs = false;
      }
      else
      {
         for(int i = 0; i < numVar; i++)
         {
            string tmp;
            pf >> tmp;
            if( pf.fail() )
            {
               cerr << "Error reading variable for index " << i << endl;
               rs = false;
         break;
            }
            var.push_back( tmp );
         }
      }
   }

   return rs;
}


bool processINT()
{
   cout << endl << "Processing INT section..." << endl;

   bool rs = false;

   string section;

   pf >> section;

   if( section != "INT" )
   {
      cerr << "INT expected.   Read instead " << section << endl;
   }
   else
   {
      auto numInt = 0;
   
      pf >> numInt; // number of variables
   
      if( pf.fail() || numInt < 0 )
      {
         cerr << "Invalid number after INT" << endl;
      }
      else
      {
         isInt.resize(var.size());
   
         for( auto it = isInt.begin(); it != isInt.end(); ++it )
         {
            *it = false;
         }
      
         if( numInt > 0) {
            // cout << "Integer variables:" << endl;
            for( int i = 0; i < numInt; ++i )
            {
               int idx;

               pf >> idx;
               if( pf.fail() )
               {
                  cerr << "Error reading integer index " << i << endl;
            goto TERMINATE;
               }
               isInt[idx] = true;
            }
         }
         rs = true;
      }
   }

TERMINATE:
   return rs;
}


bool processCON()
{
   cout << endl << "Processing CON section..." << endl;

   bool rs = false;

   string section;

   pf >> section;

   if( section != "CON" )
   {
      cerr << "CON expected.   Read instead " << section << endl;
   }
   else
   {
      pf >> numCon >> numBnd; // we don't really do anything with numBnd

      if( pf.fail() || numCon < 0 || numBnd < 0 )
      {
         cerr << "Invalid number(s) after CON" << endl;
      }
      else
      {
      
         for( int i = 0; i < numCon; i++ )
         {
            SVectorGMP coef;
      
            string label;
            int sense;
            mpq_class rhs;
      
            rs = readConstraint( label, sense, rhs, coef );
      
            if (!rs) break;
      
            constraint.push_back( Constraint( label, sense, rhs, coef, false, emptyList) );
         }
      }
   }

   return rs;
}


bool processOBJ()
{
   cout << endl << "Processing OBJ section..." << endl;

   bool rs = false;

   string section;

   pf >> section;

   if( section != "OBJ" )
   {
      cerr << "OBJ expected.   Read instead " << section << endl;
   }
   else
   {
      SVectorGMP coef;
      string objsense;
      mpq_class rhs;

      pf >> objsense;

      if( objsense == "min" )
      {
          isMin = true;
      }
      else if( objsense == "max" )
      {
          isMin = false;
      }
      else
      {
          cerr << "Invalid objective sense: " << objsense << endl;
          goto TERMINATE;
      }

      rs = readConstraintCoefs( objCoef );

      if (!rs)
      {
         cerr << "Failed to read objective coefficients" << endl;
      }
   }

TERMINATE:
   return rs;
}


bool processSOL()
{
   cout << endl << "Processing SOL section..." << endl;

   bool rs = false;
   mpq_class val;

   string section, label;

   pf >> section;

   if( section != "SOL" )
   {
      cerr << "SOL expected.   Read instead " << section << endl;
      return rs;
   } 
   
   pf >> numSol;

   if( pf.fail() )
   {
      cerr << "Failed to read number after SOL" << endl;
   }
   else if( numSol < 0 )
   {
      cerr << "Invalid number after SOL: " << numSol << endl;
   }
   else
   {
      auto satisfies = [] (Constraint &con, SVectorGMP &x)
      {
         bool rstat = false;

         mpq_class prod = scalarProd( con.coefSVec(), x );

         if( con.getSense() < 0 )
         {
            rstat = ( prod <= con.getRhs() );
         }
         else if( con.getSense() > 0 )
         {
            rstat = ( prod >= con.getRhs() );
         }
         else
         {
            rstat = ( con.getRhs() == prod );
         }

         return rstat;
      };

      SVectorGMP solSp;
      vector<mpq_class> sol( numVar );

      for( int i = 0; i < numSol; ++i)
      {
         pf >> label;
         cout << "checking solution " << label << endl;
   
         if( !readConstraintCoefs( solSp ) )
         {
            cerr << "Failed to read solution." << endl;
            goto TERMINATE;
         }
         else
         { 
            for( int j = 0; j < numVar; ++j)
            {
               sol[j] = 0;
            }

            // check integrality constraints 
            for( auto it = solSp.begin(); it != solSp.end(); ++it)
            {
               if( isInt[it->first] && !isInteger( it->second ) )
               {
                  cerr << "Noninteger value for integer variable "
                       << it->first << endl;
                  goto TERMINATE;
               }
               sol[it->first] = it->second;
            }

            for( int j = 0; j < numCon; ++j)
            {
               if ( !satisfies( constraint[i], solSp ) )
               {
                  cerr << "Constraint " << i << " not satisfied." << endl;
                  goto TERMINATE;
               }
            }
         }
         
         val = scalarProd( objCoef , solSp );

         cout << "   objval = " << val << endl;

         // update best obj fn val
         if( i )
         {
            if( isMin && val < bestObjVal )
            {
               bestObjVal = val;
            }
            else if( !isMin && val > bestObjVal )
            {
               bestObjVal = val;
            }
         }
         else
         {
            bestObjVal = val;
         }
      }

      if( numSol )
      {
         cout << "Best objval: " << bestObjVal << endl;

         if( isMin && checkUpper && bestObjVal > upperBound )
         {
            cerr << "Upper bound violated." << endl;
            goto TERMINATE;
         }
         else if( !isMin && checkLower && bestObjVal < lowerBound )
         {
            cerr << "Lower bound violated." << endl;
            goto TERMINATE;
         }
      }

      rs = true;
   }

TERMINATE:
   return rs;
}


bool processDER()
{

   cout << endl << "Processing DER section..." << endl;

   bool rs = false;

   string section;

   pf >> section;

   if( section != "DER")
   {
      cerr << "DER expected.   Read instead " << section << endl;
      return false;
   }
   
   pf >> numDer;

   cout << "Number of constraints to be derived: " << numDer << endl << endl;

   for( int i = 0; i < numDer; ++i )
   {

      string label;
      int sense;
      mpq_class rhs;

      SVectorGMP coef;

      if( !readConstraint( label, sense, rhs, coef ))
      {
         return false;
      }

      // obtain derivation method and info
      string bracket, kind;
      int refIdx;

      pf >> bracket >> kind;

      if( bracket != "{" )
      {
         cerr << "Expecting { but read instead " << bracket << endl;
         return false;
      }

      DerType derType = DerType::UNKNOWN;

      if( kind == "asm" )
         derType = DerType::ASM;
      else if( kind == "lin" )
         derType = DerType::LIN;
      else if( kind == "rnd" )
         derType = DerType::RND;
      else if( kind == "uns" )
         derType = DerType::UNS;

      // the constraint to be derived
      Constraint toDer(label, sense, rhs, coef, (derType == DerType::ASM), emptyList); 

      cout << numCon + i << " - deriving..." << label << endl;

      SVectorBool asmList;

      int newConIdx = constraint.size();

      switch (derType)
      {
         case DerType::ASM:
            asmList[ newConIdx ] = true;
            pf >> bracket;

            if( bracket != "}" )
            {
               cerr << "Expecting } but read instead " << bracket << endl;
               return false;
            }
            break;
         case DerType::LIN:
         case DerType::RND:
            {
               SVectorGMP coefDer;
               mpq_class rhsDer;
               int senseDer;

              if( !readLinComb( senseDer, rhsDer, coefDer, newConIdx, asmList) )
                 break;

               pf >> bracket;

               if( bracket != "}" )
               {
                  cerr << "Expecting } but read instead " << bracket << endl;
                  return false;
               }
            
               Constraint derived("", senseDer, rhsDer, coefDer, toDer.isAsm(), 
                                           toDer.getAsmList() );


               if (derType == DerType::RND)   // round the coefficients
                  if ( !derived.round() )
                     return false;

   
               // check the actually derived constraint against the given
               if( !derived.dominates( toDer ) )
               {
                  cout << "Failed to derive constraint " << label << endl;
                  toDer.print();
   
                  cout << "Derived instead " << endl;
                  derived.print();
                  return false;
               }
            }
            break;
         case DerType::UNS:
            {
               int con1, asm1, con2, asm2;

               pf >> con1 >> asm1 >> con2 >> asm2;

               if( pf.fail() )
               {
                  cerr << "Error reading con1 asm1 con2 asm2" << endl;
                  return false;
               }

               if( (con1 < 0) || (con1 >= newConIdx) )
               {
                  cerr << "con1 out of bounds: " << con1 << endl;
                  return false;
               }

               if( (con2 < 0) || (con2 >= newConIdx))
               {
                  cerr << "con2 out of bounds: " << con2 << endl;
                  return false;
               }

               if( !canUnsplit( toDer, con1, asm1, con2, asm2, asmList ) )
               {
                  cerr << label << ": unsplit failed" << endl;
                  return false;
               }

               pf >> bracket;
               if( bracket != "}" )
               {
                  cerr << "Expecting } but read instead " << bracket << endl;
                  return false;
               }
            }
            break;
         default:
            cout << label << ": unknown derivation type " << kind << endl;
            return false;
            break;
      }

      // set the list of assumptions
      toDer.setAsmList( asmList );

      pf >> refIdx;
      toDer.setMaxRefIdx( refIdx );
      constraint.push_back( toDer );

      if( i < numDer - 1 ) // Never trash last constraint
         if( (refIdx >= 0) && (refIdx < int(constraint.size())) ) 
         {
             constraint.back().trash();
         }

#ifndef NDEBUG
      toDer.print();
#endif
   }

   cout << endl;

   auto asmList = constraint.back().getAsmList();

   if( asmList != emptyList )
   {
      cout << "Final derived constraint undischarged assumptions:" << endl;
      for( auto it = asmList.begin(); it != asmList.end(); ++it )
      {
         auto idx = it->first;

         cout << idx << ": " << constraint[idx].label() << endl;
      }
   }
   else
   {
      if( rtpType == RtpType::INFEAS)
      {
         if( constraint.back().isFalsehood() )
         {
            cout << "Infeasibility verified." << endl;
            rs = true;
         }
         else
            cout << "Failed to verify infeasibility." << endl;
      }
      else if( (isMin && checkLower) || (!isMin && checkUpper) )
      {
         if( rtp.isTautology() )
         {
            cout << "RTP is a tautology." << endl;
            rs = true;
         }
         else if( !constraint.back().dominates( rtp ) )
         {
            if (isMin)
               cout << "Failed to derive lower bound." << endl;
            else
               cout << "Failed to derive upper bound." << endl;
         }
         else
         {
            if (numSol) {
               cout << "Best objval over all solutions: " << bestObjVal << endl;
            }

            cout << "Successfully verified optimal value range "
                   << (lowerStr == "-inf" ? "(" : "[")
                   << lowerStr << ", " << upperStr
                   << (upperStr == "inf" ? ")" : "]")
                   << "." << endl;

            rs = true;
         }
      }
   }

   return rs;
} // processDER


inline mpq_class floor( const mpq_class &q )
{
   mpz_class z = q.get_num() / q.get_den();
   return mpq_class(z);
}


inline mpq_class ceil( const mpq_class &q )
{
   mpz_class z = (q.get_num() + q.get_den()-1) / q.get_den();
   return mpq_class(z);
}


bool isInteger( const mpq_class &q )
{
   return ( q == floor(q) );
}


bool readLinComb( int &sense, mpq_class &rhs, SVectorGMP &coef, int currConIdx,
                  SVectorBool &asmList )
{
   bool rs = true;

   SVectorGMP mult;

   if( !readMultipliers(sense, mult ) )
   {
      rs = false;
   }
   else
   {
      rhs = 0;
      coef.clear();
      asmList.clear();
      mpq_class t;
   
      for( auto it = mult.begin(); it != mult.end(); ++it )
      {
         auto idx = it->first;
         auto a = it->second;
   
         auto myAsmList = constraint[idx].getAsmList();
   
         for( auto it2 = myAsmList.begin(); it2 != myAsmList.end(); ++it2 )
            asmList[it2->first] = true;
   
         Constraint &con = constraint[idx];

         if ( con.isTrashed() ) 
         {
            cerr << "Accessing trashed constraint: " << con.label() << endl;
            rs = false;
         }
         else
         {
            SVectorGMP c = constraint[idx].coefSVec();

            for( auto itr = c.begin(); itr != c.end(); ++itr )
               coef[ itr->first ] += a * itr->second;
   
            rhs += a * constraint[idx].getRhs();
   
            if ( (constraint[idx].getMaxRefIdx() <= currConIdx) &&
                (constraint[idx].getMaxRefIdx() >= 0 ) )
               constraint[idx].trash();
         }
      }
   }

   return rs;
}


bool readMultipliers( int &sense, SVectorGMP &mult )
{

   int k;
   bool rs = true;

   sense = 0;
   mult.clear();

   pf >> k;

   for( auto j = 0; j < k; ++j )
   {
      mpq_class a;
      int idx;

      pf >> idx >> a;

      if( a == 0 ) continue; // ignore 0 multiplier

      mult[idx] = a;

      if( sense == 0 )
      {
         sense = constraint[idx].getSense() * sgn(a);
      }
      else
      {
         int tmp = constraint[idx].getSense() * sgn(a);
         if( tmp != 0 && sense != tmp)
         {
            cerr << "Coefficient has wrong sign for index " << idx << endl;
            rs = false;
            goto TERMINATE;
         }
      }
   }

TERMINATE:
   return rs;
}


bool readConstraintCoefs( SVectorGMP &v )
{
   auto rs = false;
   int k = 0;
   string tmp;

   v.clear();
   
   pf >> tmp;

   if( tmp == "OBJ" )
   {
      v = objCoef;
      rs = true;
   }
   else
   {
      k = atoi(tmp.c_str());

      if( pf.fail() )
      {
         cerr << "Error reading number of elements " << endl;
         goto TERMINATE;
      }
      else
      {
         for(int j = 0; j < k; j++)
         {
            int idx;
            mpq_class a;
         
            pf >> idx >> a;
            if( pf.fail() )
            {
               cerr << "Error reading integer-rational pair " << endl;
               goto TERMINATE;
            }
            else if( idx < 0 || idx >= numVar )
            {
               cerr << "Index out of bounds: " << idx << endl;
               goto TERMINATE;
            }
            v[idx] = a;
         }
         rs = true;
      }
   }

   v.compactify();

TERMINATE:
   return rs;
}


bool readConstraint( string &label, int &sense, mpq_class &rhs, 
                     SVectorGMP &coef )
{

   auto rs = false;
   char senseChar;

   pf >> label >> senseChar;

   if (!pf.fail())
   {
      if( senseChar == 'E' )
         sense = 0;
      else if( senseChar == 'L')
         sense = -1;
      else if( senseChar == 'G')
         sense = 1;
      else
      {
        cerr << "Unknown sense for " << label << ": " << senseChar << endl;
        goto TERMINATE;
      }
      
      pf >> rhs;

      if( !pf.fail() )
         rs = readConstraintCoefs( coef );

      if( !rs ) cerr << label <<   ": Error reading constraint " << endl;
   }

TERMINATE:
   return rs;
}


// con1 and con2 must be inequalities for an integer disjunction
// e.g. mx <= d and mx >= d+1 such that the variables indexed by
// the support of m are integers.   The function checks this.
bool canUnsplit( Constraint &toDer, const int con1, const int a1, 
                 const int con2, const int a2, SVectorBool &asmList )
{

   bool rs = false;

   Constraint &c1 = constraint[con1];
   Constraint &c2 = constraint[con2];

   if (c1.isTrashed())
   {
      cerr << "unsplitting trashed constraint: " << c1.label() << endl; 
      goto TERMINATE;
   } 
   else if (c2.isTrashed())
   {
      cerr << "unsplitting trashed constraint: " << c2.label() << endl; 
      goto TERMINATE;
   }

   if( c1.dominates( toDer) && c2.dominates( toDer ) )
   {
      SVectorGMP asm1Coef, asm2Coef;
      mpq_class asm1Rhs, asm2Rhs;
   
      SVectorBool asm1 = c1.getAsmList();
      SVectorBool asm2 = c2.getAsmList();

      // remove the indices involved in unsplitting
      asm1.erase( a1 );
      asm2.erase( a2 );

      asmList.clear();
      asmList = asm1;

#ifndef NDEBUG
      cout << "asm1: ";
      for(auto it = asm1.begin(); it != asm1.end(); ++it) {
         cout << it->first << " ";
      }
      cout << endl;

      cout << "asm2: ";
      for(auto it = asm2.begin(); it != asm2.end(); ++it) {
         cout << it->first << " ";
      }
      cout << endl;
#endif

      for(auto it = asm2.begin(); it != asm2.end(); ++it) {
         if (asmList.find(it->first) != asmList.end()) {
            asmList[it->first] = true;
         }
      }

      c1 = constraint[a1];
      c2 = constraint[a2];

      if (c1.isTrashed())
      {
         cerr << "accessing trashed constraint: " << c1.label() << endl; 
         goto TERMINATE;
      } 
      else if (c2.isTrashed())
      {
         cerr << "accessing trashed constraint: " << c2.label() << endl; 
         goto TERMINATE;
      }


      // the constraints must have opposite senses
      if( -1 != c1.getSense() * c2.getSense() )
      {
         cerr << "canUnsplit: Failed sense requirement for assumptions" << endl;
         cerr << "c1 sense:: " << c1.getSense() << endl;
         cerr << "c2 sense:: " << c2.getSense() << endl;
         goto TERMINATE;
      }
      else
      {
         // check if disjunction gives a tautology with respect to the variable
         // integrality requirements

         bool stat = true;
         if( c1.getSense() < 0 )
            stat = ((c1.getRhs() + 1) == c2.getRhs());
         else // must be > 0
            stat = (c1.getRhs() == (c2.getRhs() + 1));

         if( !stat ) {
            cerr << c1.label() << " and " << c2.label()
                 << " do not form a tautology" << endl;
            goto TERMINATE;
         };

         if( c1.coefSVec() != c2.coefSVec() )
         {
            cerr << "canUnsplit: coefs of asm constraints differ" << endl;
            goto TERMINATE;
         }
         else
         {
            SVectorGMP c = c1.coefSVec();
            
            for( auto it = c.begin(); it != c.end(); ++it)
            {
               if( !isInt[it->first] )
               {
                  cerr << "canUnsplit: noninteger variable index " << it->first
                         << endl;
                  goto TERMINATE;
               }
               else if( !isInteger(it->second) )
               {
                  cerr << "canUnsplit: noninteger coefficient for index " 
                         << it->first << endl;
                  goto TERMINATE;
               }
            }
         }
         rs = true;
      }
   }

TERMINATE:
   return rs;
}


bool SVectorGMP::operator!=( SVectorGMP &other )
{
   bool rs = false;

   SVectorGMP &cf1 = *this;
   SVectorGMP &cf2 = other;

   // get rid of all zero entries
   cf1.compactify();
   cf2.compactify();

   if (cf1.size() != cf2.size())
   {
      rs = true;
   }
   else
   {
      for(auto it1 = cf1.begin(); it1 != cf1.end(); ++it1)
      {
         auto it2 = cf2.find(it1->first);
         if( it2 == cf2.end() || it2->second != it1->second )
         {
            rs = true;
            break;
         }
      }
   }

   return rs;
}


// use non-sparse vector for v to reduce lookup time
mpq_class scalarProd( const SVectorGMP &u, const SVectorGMP &v )
{
   mpq_class prod = 0;

   for( auto it = u.begin(); it != u.end(); ++it )
   {
      auto it2 = v.find( it->first );
      if( it2 != v.end() )
         prod += it->second * it2->second;
   }

   return prod;
}

bool Constraint::round()
{
   bool rs = true;

   for( auto it = _coef.begin(); it != _coef.end(); ++it )
   {
      auto j = it->first;
      auto a = it->second;

      if( isInt[j] )
      {   // needs to be an integer variable
         if( !isInteger( a ) )
         {
            cerr << "Coefficient of integer variable with index "
                 << j << " is not an integer" << endl;
            rs = false;
            goto TERMINATE;
         }
      }
   }

   if( getSense() < 0) // round down
      _rhs = floor( _rhs );
   else if (getSense() > 0) // round up
      _rhs = ceil( _rhs );

TERMINATE:
   return rs;
}


bool Constraint::_isFalsehood()
{
   bool rs = false;

   if( _coef.size() == 0 )
   { 
      if(   ((getSense() <= 0) && (_rhs < 0)) 
          || ((getSense() >= 0) && (_rhs > 0)) )
         rs = true;
   }

   return rs;
}


bool Constraint::dominates( Constraint &other )
{ 
   bool rs = false;

   if( this->_isFalsehood() )
   {
      rs = true;
   }
   else if (this->_coef == other._coef) 
   {
      if(   (other.getSense() > 0 && this->getSense() >= 0 && 
             this->_rhs >= other._rhs)
         || (other.getSense() < 0 && this->getSense() <= 0 &&
             this->_rhs <= other._rhs)
         || (other.getSense() == 0 && this->getSense() == 0 &&
             this->_rhs == other._rhs) )
      {
         rs = true;
      }
   }

   return rs;
}

 
bool Constraint::isTautology() {
   bool rs = false;

   if (_coef.size() == 0) {
      if (    ((getSense() == 0) && (0 == _rhs))
            || ((getSense() < 0) && (_rhs >= 0 ))
            || ((getSense() > 0) && (_rhs <= 0 ))) {
         rs = true;
      }
   }
   return rs;
}


void Constraint::print() {
   bool first = true;
   mpq_class myCoef;

   int cnt = 0;

   if (_isAsm)
      cout << "Is assumption: ";

   for(auto it = _coef.begin(); it != _coef.end(); ++it)
   {
      auto idx = it->first;
      auto a = it->second;

      myCoef = abs(a);

      if( a > 0 )
      {
         if( !first ) cout << " + ";
         cout << (myCoef == 1 ? string("") : string(myCoef.get_str()) + " ") 
                << var[idx];
         ++cnt;
         first = false;
      }
      else if( a < 0 )
      {
         cout << " - " 
                << (myCoef == 1 ? string("") : string(myCoef.get_str()) + " ") 
                << var[idx];
         first = false;
         ++cnt;
      }

      if( (cnt+1) % 4 == 0 ) cout << endl;
   }

   if( first ) // coefficients are all zero
   { 
      cout << "0";
   }

   switch( _sense )
   {
      case -1: cout << " <= "; break;
      case 1: cout << " >= "; break;
      case 0: cout << " = "; break;
      default : assert(false);
   }
   cout << _rhs << endl;

#ifndef NDEBUG
   if( !_isAsm)
   {
      cout << " -- assumptions: " << endl;
      for( auto it = _asmList.begin(); it != _asmList.end(); ++it )
      {
         auto idx = it->first;
         cout << "   "<< it->first << ": " << constraint[idx].label() << endl;
      }
      cout << endl;
   }
#endif
}

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
#include <vector>
#include <set>

#define VERSION_MAJOR 1
#define VERSION_MINOR 0

using namespace std;



int main(int argc, char *argv[])
{

   int fpos;
   bool toTrim = false;
   int trimmed = 0;

   ifstream pf; // input vipr file
   ofstream tnf;  // tightened/trimmed file

   vector<int> maxConIdx;
   vector<int> newConIdx;
   vector<set<int>> neededBy;
   vector<set<int>> needs;


   int rs = -1;
   int farg = 1;

   if( (argc == 1) || (argc > 3) )
   {
      cerr << "Usage: " << argv[0] << " [-t] filename\n";
      cerr << "  Specify option -t to trim unused derived constraints" << endl;
      return rs;
   }
   else if( argc == 3 )
   {
      if( string(argv[1]) == "-t" ) 
      {
         toTrim = true;
         farg = 2;
      }
      else if( string(argv[2]) == "-t" )
      { 
         toTrim = true;
         farg = 1;
      }
      else
      {
         cerr << "Usage: " << argv[0] << " [-t] filename\n";
         cerr << "  Specify option -t to trim unused derived constraints" 
              << endl;
         return rs;
      }  
   }


   pf.open( argv[farg] );

   if( pf.fail() )
   {
      cerr << "Failed to open file " << argv[farg] << endl;
      return rs;
   }

   string tnFname = string(argv[farg]) 
                       + string(toTrim ? ".trimmed" : ".tightened");

   tnf.open( tnFname.c_str());

   if ( tnf.fail() )
   {
      cerr << "Failed to open file " << tnFname << endl;
   }
   else
   {

      auto _checkVersion = [](string ver)
      {
         bool rstat = false;
         int major, minor;
      
         size_t pos = ver.find( "." );
      
         major = atoi( ver.substr( 0, pos ).c_str() );
         minor = atoi( ver.substr( pos+1, ver.length()-pos ).c_str() );
      
         if ( (major == VERSION_MAJOR) && (minor <= VERSION_MINOR ) )
         {
            rstat = true;
         }
         else
         {
            cerr << "Version " << ver << " unsupported" << endl;
         }
      
         return rstat;
      };

      auto _processSparseVec = [ &pf, &tnf, &newConIdx ](bool toPrint, bool useNewIdx)
      {
         bool rval = true;

         string input, val;
         int k, index;

         pf >> input;
         if (input == "OBJ")
         {
            if (toPrint) tnf << " OBJ ";
         }
         else
         {
            k = atoi(input.c_str());                

            if (toPrint) tnf << " " << k;

            for( int i = 0; i < k; ++i )
            {
               pf >> index >> val;
               if( pf.fail() )
               {
                  cerr << "Failed reading coefficient " << i << endl;
                  rval = false;
            break;
               }
               else if (toPrint)
               {
       
                  if ( (i+1) % 10 == 0) tnf << endl;

                  if( useNewIdx ) index = newConIdx[index];

                  tnf << "  " << index << " " << val;
               }
            }
         }
         return rval;
      };

      auto _updateMaxConIdx = [ &pf, &tnf, &maxConIdx, &needs, &neededBy ](int curConIdx)
      {
         bool rval = true;

         string input, val;
         int k, index;

         pf >> input;
         k = atoi(input.c_str());                

         for( int i = 0; i < k; ++i )
         {
            pf >> index >> val;
            if( pf.fail() )
            {
               cerr << "Failed reading coefficient " << i << endl;
               rval = false;
         break;
            }
            else 
            {
               maxConIdx[ index ] = curConIdx;
               needs[curConIdx].insert( index );
               neededBy[ index ].insert( curConIdx );
            }
         }
         return rval;
      };

      string section, tmp, label;
      char sense;
      int num, idx, numCon, numSol, numDer, numBnd;
      int numAll; // total number of contraints and derived constraints
      int con1, asm1, con2, asm2; // for reading unsplitting indices
      bool stat = false;

      // Eat up comment lines, if any, until hitting VER
      for(;;)
      {
         pf >> section;

         if( pf.fail() ) goto TERMINATE;

         if (section == "VER")
         {
            pf >> tmp;
            if( _checkVersion( tmp ))
            {
               tnf << "VER " << tmp << endl;
      break;
            }
            else
            {
               goto TERMINATE;
            }
         }
         else if (section == "%")
         {
            getline( pf, tmp );
            tnf << "% " << tmp << endl;
         }
         else
         {
            cerr << endl << "\% or VER expected. Read instead "
                   << section << endl;
            goto TERMINATE;
         }
      }


      pf >> section;
      if( section != "VAR" )
      {
         cerr << "VAR expected. Read instead: " << section << endl;
         goto TERMINATE;
      }

      pf >> num;

      if( pf.fail() ) goto TERMINATE;

      tnf << "VAR " << num << endl;

      for( int i = 0; i < num; ++i )
      {
         pf >> tmp; // read variable name
         tnf << tmp << endl;
      }


      pf >> section;
      if( section != "INT" )
      {
         cerr << "INT expected. Read instead: " << section << endl;
         goto TERMINATE;
      }

      pf >> num;

      if( pf.fail() )
      {
         cerr << "Failed to read number after INT" << endl;
         goto TERMINATE;
      }

      tnf << "INT " << num << endl;

      for( int i = 0; i < num; ++i )
      {
         pf >> idx;
         if( pf.fail() ) goto TERMINATE;

         if ((i+1) % 10 == 0) tnf << endl;
         tnf << " " << idx;
      }
      tnf << endl;


      pf >> section;
      if( section != "OBJ" )
      {
         cerr << "OBJ expected. Read instead: " << section << endl;
         goto TERMINATE;
      }

      pf >> tmp;
      if( (tmp != "min") && ( tmp != "max" ) )
      {
         cerr << "Unrecognized string after OBJ: " << tmp << endl;
         goto TERMINATE;
      }

      tnf << "OBJ " << tmp << endl;

      stat = _processSparseVec( true, false );
      tnf << endl;

      if( !stat ) goto TERMINATE;


      pf >> section;
      if( section != "CON" )
      {
         cerr << "CON expected. Read instead: " << section << endl;
         goto TERMINATE;
      }

      pf >> numCon >> numBnd;

      tnf << "CON " << numCon << " " << numBnd << endl;

      for( int i = 0; i < numCon; ++i )
      {
         pf >> label >> sense >> tmp;

         tnf << label << " " << sense << " " << tmp;

         stat = _processSparseVec( true, false );
         tnf << endl;
         if( !stat ) break;

      }


      pf >> section;
      if( section != "RTP" )
      {
         cerr << "RTP expected. Read instead: " << section << endl;
         goto TERMINATE;
      }

      pf >> tmp;

      tnf << "RTP " << tmp;
      if( tmp == "range" )
      {
         pf >> tmp; // lower bound
         tnf << " " << tmp;

         pf >> tmp; // upper bound
         tnf << " " << tmp;
      }
      else if ( tmp != "infeas" )
      {
         cerr << "Unrecognized string after RTP: " << tmp << endl;
         goto TERMINATE;
      }

      tnf << endl;

      pf >> section;
      if( section != "SOL" )
      {
         cerr << "SOL expected. Read instead: " << section << endl;
         goto TERMINATE;
      }

      pf >> numSol;

      tnf << "SOL " << numSol << endl;

      if (numSol)
      {
         for( int i = 0; i < numSol; ++i ) {
            pf >> label;
            tnf << label;
            stat = _processSparseVec( true, false );
            tnf << endl;
         }
      }


      pf >> section;
      if( section != "DER" )
      {
         cerr << "DER expected. Read instead: " << section << endl;
         goto TERMINATE;
      }

      pf >> numDer;

      numAll = numCon + numDer;
      maxConIdx.resize( numAll );
      needs.resize( numAll );
      neededBy.resize( numAll );

      if( toTrim ) newConIdx.resize( numAll );

      for( auto& it : maxConIdx )
         it = 0;
      maxConIdx.back() = numAll - 1; // the last derived constraint is NEVER redundant
      neededBy.back().insert( numAll - 1 );

      fpos = pf.tellg();

      // Go through derived constraints and set max constraint indices
      for( int i = 0; i < numDer; ++i )
      {
         int conIdx = i + numCon;

         pf >> label >> sense >> tmp;

         stat = _processSparseVec( false, false );
         if( !stat ) 
         {
            cerr << "Error processing " << label << endl;
            goto TERMINATE;
         }

         pf >> tmp;
         if( tmp != "{" )
         {
            cerr << "'{' expected.   Reading instead: " << tmp << " in " 
                 << label << endl;
            goto TERMINATE;
         }
         else
         {
            pf >> tmp;
            if( tmp == "asm" )
            {
               pf >> tmp;
               if( tmp != "}")
               {
                  cerr << "'}' expected. Read instead: " << tmp << endl;
                  goto TERMINATE;
               }
            }
            else if( tmp == "lin" )
            {
               stat = _updateMaxConIdx( conIdx );
               if( stat )
               {
                  pf >> tmp;
                  if( tmp != "}")
                  {
                     cerr << "'}' expected. Read instead: " << tmp << endl;
                     goto TERMINATE;
                  }
               }
            }
            else if( tmp == "rnd" )
            {
               stat = _updateMaxConIdx( conIdx );
               if( stat )
               {
                  pf >> tmp;
                  if( tmp != "}")
                  {
                     cerr << "'}' expected. Read instead: " << tmp << endl;
                     goto TERMINATE;
                  }
               }
            }
            else if( tmp == "uns" )
            {
               pf >> con1 >> asm1 >> con2 >> asm2;
               if( pf.fail() )
               {
                  goto TERMINATE;
               }
               else
               {
                  maxConIdx[con1] = conIdx;
                  maxConIdx[con2] = conIdx;
                  maxConIdx[asm1] = conIdx;
                  maxConIdx[asm2] = conIdx;
                  needs[conIdx].insert( con1 );
                  needs[conIdx].insert( con2 );
                  needs[conIdx].insert( asm1 );
                  needs[conIdx].insert( asm2 );
                  neededBy[ con1 ].insert( conIdx );
                  neededBy[ con2 ].insert( conIdx );
                  neededBy[ asm1 ].insert( conIdx );
                  neededBy[ asm2 ].insert( conIdx );

                  pf >> tmp;
                  if( tmp != "}")
                  {
                     cerr << "'}' expected. Read instead: " << tmp << endl;
                     goto TERMINATE;
                  }
               }
            }
            else
            {
               cerr << "Unrecognized reason type: " << tmp << endl;
               goto TERMINATE;
            }

         }

         pf >> idx; // read off current max con index and ignore it

      }


      if( toTrim ) 
      { 
         // process backwards to peel of neededBy indices
         for( auto i = int(needs.size()) - 1; i >= numCon; --i)
         {
            if( neededBy[i].empty() )
            {
               for( auto it : needs[i] )
               {
                   neededBy[ it ].erase( neededBy[ it ].find( i ) );
               }
            }
         }

         idx = 0;
         for( size_t i = 0; i < neededBy.size(); ++i )
         {
            newConIdx[i] = idx;
            if( (i < size_t(numCon)) || !neededBy[i].empty() ) {
                ++idx;
            }
         }
         trimmed = numCon + numDer - idx;
      }

      rs = 0;

      pf.seekg( fpos );

      tnf << "DER " << numDer - trimmed << endl;

      for( int i = 0; i < numDer; ++i )
      {
         int conIdx = i + numCon;
         bool toPrint = (!toTrim || !neededBy[ conIdx ].empty());

         pf >> label >> sense >> tmp;

         if (toPrint) tnf << label << " " << sense << " " << tmp;

         stat = _processSparseVec( toPrint, toTrim );
         if( !stat ) 
         {
            cerr << "Error processing " << label << endl;
            goto TERMINATE;
         }

         pf >> tmp;
         if ( toPrint ) tnf << " " << tmp;
         if( tmp != "{" )
         {
            cerr << "'{' expected.   Reading instead: " << tmp << " in " 
                 << label << endl;
            goto TERMINATE;
         }
         else
         {

            pf >> tmp;
            if( toPrint ) tnf << " " << tmp;
            if( tmp == "asm" )
            {
               pf >> tmp;
               if( toPrint ) tnf << " " << tmp;
               
               if( tmp != "}")
               {
                  cerr << "'}' expected. Read instead: " << tmp << endl;
                  goto TERMINATE;
               }
            }
            else if( tmp == "lin" )
            {
               stat = _processSparseVec( toPrint, toTrim );
               if( stat )
               {
                  pf >> tmp;
                  if( toPrint ) tnf << " " << tmp;
                  if( tmp != "}")
                  {
                     cerr << "'}' expected. Read instead: " << tmp << endl;
                     goto TERMINATE;
                  }
               }
            }
            else if( tmp == "rnd" )
            {
               stat = _processSparseVec( toPrint, toTrim );
               if( stat )
               {
                  pf >> tmp;
                  if( toPrint ) tnf << " " << tmp;
                  if( tmp != "}")
                  {
                     cerr << "'}' expected. Read instead: " << tmp << endl;
                     goto TERMINATE;
                  }
               }
            }
            else if( tmp == "uns" )
            {
               pf >> con1 >> asm1 >> con2 >> asm2;

               if( toTrim )
               {
                  con1 = newConIdx[con1]; asm1 = newConIdx[asm1];
                  con2 = newConIdx[con2]; asm2 = newConIdx[asm2];
               }
               if( toPrint )
               {
                  tnf << " " << con1 << " " << asm1;
                  tnf << "  " << con2 << " " << asm2;
               }
               if( pf.fail() )
               {
                  goto TERMINATE;
               }
               else
               {
                  pf >> tmp;
                  if( toPrint ) tnf << " " << tmp;
                  if( tmp != "}")
                  {
                     cerr << "'}' expected. Read instead: " << tmp << endl;
                     goto TERMINATE;
                  }
               }
            }
            else
            {
               cerr << "Unrecognized reason type: " << tmp << endl;
               goto TERMINATE;
            }

         }

         pf >> idx; // read off current max con index and ignore it
         if( toPrint )
         {
            idx = maxConIdx[ conIdx ];
            if( toTrim ) idx = newConIdx[ idx ];
            tnf << " " << idx << endl;
         }

      }

      if( toTrim )
      {
          cout << "Removed " << trimmed << " of " << numDer << " derived constraints" 
               << endl;
      }

TERMINATE:
      if( !stat ) {
         cerr << "Error encountered while processing file" << endl;
      }
      tnf.close();
   }

   pf.close();

   return rs;
}

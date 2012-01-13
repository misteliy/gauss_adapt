   +-----------------------------------------------------------------------+
   |                                                                       |
   |   MVN - Generate random multivariate normal numbers                   |
   |   Version 1 - July 2006                                               |
   |                                                             MVN 1.0   |
   |                                                (c) John S. Uebersax   |
   |                 http://ourworld.compuserve.com/homepages/jsuebersax   |
   |                                                                       |
   +-----------------------------------------------------------------------+

The file mvn.zip contains the following:

             Readme.txt    this file
                mvn.exe    executable program
            example.txt    example input file

             /benchmark    a folder with benchmark input/output files


CITATION

Please cite MVN in any publications resulting from its use.  This will (1) let
other scientists replicate your work and (2) help others who may wish to
generate multivariate normal random data.  A suggested citation format is:

     Uebersax JS.  MVN program for random multivariate normal numbers.
     2006.  Available at the Statistical Method for Rater Agreement website:
     http://ourworld.compuserve.com/homepages/jsuebersax .  Accessed:
     mmm dd, yyyy.


FEATURES AND LIMITATIONS

MVN is a simple but technically solid program for generating random
multivariate normal numbers.

  * Minimalist design/interface
  * Up to one million random numbers (actually one less than that)
  * Unlimited number of variables
  * User can specify means, standard deviations and/or correlations
  * Runs in Command Prompt window on Windows 95/98/NT/2000/XP


INSTRUCTIONS FOR USE


1.  Input File

The input file specifies the run parameters.  It has six required lines and
optional lines.


Required lines

The following example shows the required lines:

     5   Number of variables
 50000   Number of random vectors
  1234   Random number seed
     0   Means are supplied (0=no, 1=yes)
     0   Standard deviations are supplied (0=no, 1=yes)
     0   Correlations are supplied (0=no, 1=yes)

Note that each line contains a numeric and a text field.  The numeric field
occupies columns 1-6 and contains an integer that ends *exactly* in column 6;
a blank field is read the same as a 0.

The text field is used for comments and is not read by the program.


Line by Line Explanation

     Line 1  specifies the number of variables.

     Line 2  specifies the number of random vectors.  Each vector contains
             a random value for each variable.

     Line 3  supplies a seed value for the random number generator.  This 
             must be an integer > 0.  Using the same seed value for the 
             same problem will produce identical samples for different 
             runs.  Using different seed values will produce different 
             random samples.

     Line 4  0 for default means or 1 if you supply them

     Line 5  0 for default standard deviations, or 1 if you supply them

     Line 6  0 for default correlations, or 1 if you supply them

If means, standard deviations, and/or correlations are not user-supplied
then default values are used.  Default values are as follows:

                      Default
     Parameter         value
     -----------      -------
     Means               0.
     Std devs            1.
     Correlations        0.

So, for instance, the example input file above would produce 50000 random 
vectors, each with a random value for five variables, with means = 0, 
standard deviations = 1, and correlations = 0.


Optional Lines

Here you can specify the means, standard deviations, and correlations of the
random data.

    * Only supply the values which you have indicated in lines 4-6 that 
      you will supply.  That is, you can supply means but not standard 
      deviations, correlations but neither means nor standard deviations, 
      etc.

    * Supply means first, then standard deviations, the correlations.

    * Begin each set of values on a new line.  Aside from that restriction,
      the format is free-field.

    * For correlations, supply only the lower triangle of the correlation 
      matrix, (without the diagonal)


Very Important:

Be careful to press the enter key after the last number of the last line.
Otherwise the end-of-line may not be marked correctly and the data might not be
read properly.  Some editors do this automatically, but some (including Notepad)
do not.

In the following example, 1000 random vectors of 3 numbers each are produced.
The means are all 100.  The standard deviations are all 15.  The correlations are
as follows:  r(v1, v2) = .7; r(v1, v3) = .5; r(v2, v3) = .4.

     3   Number of variables
  1000   Number of random vectors
    17   Random number seed
     0   Means are supplied (0=no, 1=yes)
     0   Standard deviations are supplied (0=no, 1=yes)
     0   Correlations are supplied (0=no, 1=yes)
 100 100 100
 15  15  15
 .7
 .5  .4


2. Running The Program

MVN runs in a Command Prompt window.  The simplest way to run the program is to
use Windows Explorer to navigate to the folder where the file mvn.exe is
located.  Then click the icon for mvn.exe.  This will open a Command Prompt
window with MVN running in it.

This is the simplest, but not necessarily the best way to run the program.
Alternatively one can open a Command Prompt window first, navigate to the
folder with mvn.exe, then type:

     mvn

and press the enter key.  For general instructions about how to use
Command Prompt, which is a very helpful "power user" tool, see:

     http://ourworld.compuserve.com/homepages/jsuebersax/dos.htm

Whichever method you use, you can adjust the size and appearance of
the Command Prompt window by right-clicking its title bar and selecting
Properties.


File Names

MVN will first prompt for the names of the input and output files.  If you
just press enter the default names of Input.txt and Data.txt will be 
assumed.  To use another name, enter the name (up to 60 characters) and press 
enter.  If you use the extension .txt for these files you can open them in
the Notepad editor by clicking on their icons.

You can include a path-specification along with the file name.


For Excel Users

If you choose the extension .csv for the output file it will have a 
comma-separated values format.  This means if you click on its icon the file
should open in Excel automatically.


3. Output Format

The standard output format is 12f9.4.  However, if you specify a .csv extension
for the output file the format is 10f15.6.

If neither format is satisfactory and if you would prefer a version of MVN that
lets the user specify the output format, please email me.


4. Troubleshooting

     1.  Is the input file in the same folder as mvn.exe?  Or, if not, have
         you supplied the correct path? 

     2.  When you make the input file, make sure you press the enter key after 
         the last value of the last line. 

     3.  Can you replicate the benchmark results?

     4.  Do values end in column 6 in lines 1-6?

     5.  Correct number of means/standard deviations, if supplied?

     6.  Lower triangle correlation matrix, no diagonal, if correlations 
         supplied?

     7.  If you ran MVN by clicking the icon of the file mvn.exe, and if the
         Command Prompt window closed before you could read what it said, try 
         opening the Command Prompt window first, navigating to the folder with 
         mvn.exe, type

              mvn

         and press the enter key.  This will keep the Command Prompt window 
         open.

If you supply an improper (not positive definite) correlation matrix, the 
Cholesky decomposition will fail and you will get no results.  An (extreme) 
example of an improper correlation matrix would be r(v1, v2) = 1, r(v1, v3) 
= 1 and r(v2, v3) = -1.  Note that this correlation structure is impossible. 

If the correlation matrix you supply is improper, it will be written to
the output file you specified to help you in diagnosing the problem.

If none of the above solve your problem, feel free to email me; please include
your input file.


TECHNICAL

MVN was written in Fortran 90.

To produce random multivariate normal numbers, MVN first generates random
univariate normal numbers.  This is done using the TOMS Algorithm 712 by JL
Leva.  The full reference is:

    Leva JL.  Algorithm 712. A normal random number generator.
    ACM Transactions on Mathematical Software (TOMS), v.18 n.4,
    pp. 454-455, Dec. 1992

The algorithm uses the ratio of uniforms method of AJ Kinderman and JF
Monahan augmented with quadratic bounding curves (citation needed).

Uniform random numbers, used by this algorithm, are supplied by the default
random number function of the Fortran 90 compiler (Absoft Pro Fortran 90, v. 7).

The random multivariate normal numbers are produced by pre-multiplying a vector
of random univariate normal numbers by the Cholesky decomposition of the
correlation matrix according to the formula:

     Y = L X

where

     Y = a vector of random multivariate normal numbers

     X = a vector of random univariate normal numbers

     L = the Cholesky decomposition of the correlation matrix, stored
         in the lower triangle and main diagonal of a square matrix
         (elements in the upper triangle of the matrix are 0.)

Standard deviations are then multiplied and/or means added per the user
specifications.

Happy computing!

John S. Uebersax
jsuebersax@gmail.com

§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§

Terms and Conditions

The author can make no guarantee concerning the accuracy or correct 
working of this program. The user assumes all associated risks.  It is 
recommended that the user check the random numbers produced for 
conformance to the specified means, standard deviations and 
correlations.

MVN is free software.  

§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§»«§

History:

v 1.0 (July 2006)

  - First version

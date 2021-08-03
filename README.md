# knot-confinement-and-PH
Dataset for the paper "A statistical approach to knot confinement via PH" with B.I.Mahler

This dataset consists of several .xls files, and a single .py file. The former are divided in two kinds, called "all_knots_with_size_X" 
and "specific-X-Y". 
In the first case X is an integer between 10 and 100 (with step 10); each file contains 10.000 PL knot embeddings, decorated with:
the name of the knot (as created by Topoly), the diameter of the minimal circumscribing sphere, the knot's curvature and torsion,
the integral of the VR first Betti curve, the number of bars appearing in the VR barcode (in dimension 1), the length of the maximal
bar, the average crossing number, the writhe, the volume of the convex hull and the radius of gyration.

The same data is collected in the "specific-X-Y", but this time X is a knot type (for knots with less than 7 crossings, and including 
an "unknown" category that collects knots with minimal crossing number greater than 6), and Y is an integer between 50 and 200 with 
step 50. For each pair X,Y we examine 1000 PL knots.

Finally, the Python file "github_python_functions.py" contains the main functions we used to produce the dataset described above.

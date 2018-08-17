# Quantization-and-QM-coder

QM encoder is an adaptive binary arithmetic coding procedure which is widely used in entropy coding such as JPEG. 
It involves sub division of intervals recursively. I’ve implemented QM coding using the algorithm given, and has compared 
the performance of QM encoding on various data files.

Approach and Procedure
QM coding essentially involves loading the given transition table and updating and retrieving the current state, Qe (probability of LPS) 
and thereby using the algorithm to update the values of A and C intervals. In the process of executing QM coding, 
I’ve created and used 3 classes along with headers (one to load and retrieve appropriate values from QM table, 
another to handle file input output operation, another class which is called if mapping is required). As a part of mapping the data, 
I’ve used bit plane mapping technique. It is a process of putting together all MSBs, all LSBs and other bits in the same manner. 
I’ve defined a user input in main function which asks the user to whether apply mapping or not and acts accordingly. 
After successfully opening the input and output files, I’ve initialized the values corresponding to state transition table 
as per the requirements. Then I’ve read the data and have checked whether it is a MPS or LPS and performed the respective 
operations as and when it is necessary. Whenever the required condition is met (during renormalization and carry) 
I’ve written the corresponding bits to the output file and has reported the file size after the process is done. There after 
I’ve computed the compression ratio.

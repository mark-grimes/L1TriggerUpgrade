/*
/  README: Describes format of the L1Menu.txt file and what is required to
/          add a new trigger.
/
*/

1) L1Menu*.txt is a text file that describes the possible L1 Algorithms
   that **can be** used in L1Menu2015.C. The "can be" is highlighted because
   as described below, algorithms can be deactivated by a flag in the
   text file.  This file is read in by L1Menu2015.C (InitL1Menu Method) prior 
   to running over a sample.
	
2) Format of L1Menu.txt File:

    Column 1: (string) Name of algorithm. Should start with "L1_" and be less
                       than 20 characters.
					
    Column 2: (int)    Trigger bit for algorithm.  Not critical but keeps
                       formating inline with older code.
							  
    Column 3: (int)    Prescale value.  Nominally would be "1". IF SET TO
                       ZERO, THAT DISABLES THE ALGORITHM AND IT IS NOT USED
                       IN THE MENU EVALUATION.
							  
    Column 4: (float)  Primary Et/Pt threshold
    Column 5: (float)  Secondary Et/Pt threshold (if used)
    Column 6: (float)  Third Et/Pt threshold (if used)
    Column 7: (float)  Forth Et/Pt threshold (if used)
	 
    Column 8: (float)  Eta cut. Which object(s) this is applied to depends
                       on the implementation in L1Menu2015.C.  Care must be
                       taken when cutting on muons vs calorimeter objects. 
                       Muon cuts tend to be in true eta value, while cal cuts
                       tend to be on eta **bin**.
							  
    Column 9: (int)    Minimum Muon Quality cut.
	 
    Column 10: (float) Allocated bandwidth.  This value is NOT used by 
                       L1Menu2015.C.  It is used by EvaluateL1Menu.C to 
                       determine new threshold values.
							  
    Column 11: (int)   Flag allowing bandwidth to be scaled up or down. Value
                       =1 allows scaling. Value = 0 fixes the bandwidth. This
                       value is NOT used by L1Menu2015.C.  It is used by 
                       EvaluateL1Menu.C to determine new threshold values.
							  	 
    Column 12: (int)   Flag allowing threshold to be locked. Value
                       =1 locks threshold. Value = 0 will allow the
                       theshold to be changed. This value is NOT used 
                       by L1Menu2015.C.  It is used by EvaluateL1Menu.C
                       to determine new threshold values. In this case, 
                       if the threshold is locked, EvaluateL1Menu.C will
                       not be allowed to change it.
							  
3) Corresponding Code in L1Menu2015.C:

  a)For each algorithm that is defined in L1Menu.txt, there must be a 
    corresponding call in the EvalMenu method of L1Menu2015.C.  In the
    method, a call to "InsertInMenu" must exist.  The name of the algorithm
    must match what is in L1Menu.txt. Note, the "InsertInMenu" method 
    depends on a boolean function that actually implements the trigger
    algorithm.                              
	

  b)If the EvaluateL1Menu.C is going to used to adjust thresholds to fill
    out a menu, it depends on a "rate vs threshold" plot.  These plots are
    made in L1Menu2015.C  in the EvalThresh method.  If a new algorithm
    is added, there must be a corresponding "rate vs threshold" plot added
    to the EvalThreshold method.  The naming convention of histograms in
    EvalThreshold should be obeyed (ie. "h_(substring alg name)_byThreshold".

#!/usr/bin/env python

# ##################################################################################################
#  Disclaimer                                                                                      #
#  This file is a python3 translation of AutoDockTools (v.1.5.7)                                   #
#  Modifications made by Valdes-Tresanco MS (https://github.com/Valdes-Tresanco-MS)                #
#  Tested by Valdes-Tresanco-MS and Valdes-Tresanco ME                                             #
#  There is no guarantee that it works like the original distribution,                             #
#  but feel free to tell us if you get any difference to correct the code.                         #
#                                                                                                  #
#  Please use this cite the original reference.                                                    #
#  If you think my work helps you, just keep this note intact on your program.                     #
#                                                                                                  #
#  Modification date: 2/5/20 19:51                                                                 #
#                                                                                                  #
# ##################################################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/calc_energy_breakdown_from_results.py,v 1.2.2.1 2015/08/26 22:44:52 sanner Exp $
#

wt_dict = {}
wt_dict['gauss 1'] = -0.035579
wt_dict['gauss 2'] = -0.005156
wt_dict['repulsion'] = 0.840245
wt_dict['hydrophobic'] = -0.035069
wt_dict['Hydrogen'] = -0.587439
wt_dict['rot'] = 0.05846
wt_dict['num_tors_div'] = 1.923

verbose = False
vina_score_only_result_file = None
vina_num_torsions = 0
intramolecular_energy = 0

if __name__ == '__main__':
    import sys
    import getopt


    def usage():
        print("Usage: calc_energy_breakdown_from_vina_output.py -f vina_score_only_result.txt")
        print()
        print("    Description of command...")
        print("         -f     vina_score_only_result.txt ")
        print("         -n     number of torsions allowed by vina ")
        print("                (vina disallows bonds rotating hydroxyls) ")
        print("         -i     intramolecular energy ")
        print("    Optional parameters:")
        print("        [-v]    verbose output")
        print("        [-o energy_results ('.txt')] (default output filename is energy_results.txt)")

    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'vf:n:i:oh')
        print("in calc_energy_breakdown_from_vina_output.py")      
    except getopt.GetoptError as msg:
        print('calc_energy_breakdown_from_vina_output.py: %s' %msg)
        usage()
        #sys.exit()

    #vina_score_only_result_file = "result_file_vina_score_only.txt"
    for o, a in opt_list:
        #print "o=", o, " a=", a
        if o in ('v', '--v'):
            verbose = True
            if verbose: print('verbose=', verbose)
        if o in ('-f', '--f'):
            vina_score_only_result_file = a
            if verbose: print('set vina_score_only_result_file to ', a)
        if o in ('-n', '--n'):
            vina_num_torsions = int(a)
            if verbose: print('set vina_num_torsions to ', a)
        if o in ('-i', '--i'):
            intramolecular_energy = float(a)
        if o in ('o', '--o'):
            output_filename = a
            if verbose: print('set energy_results file to ', a)

    if not vina_score_only_result_file:
        print('calc_energy_breakdown_from_vina_output.py: vina_score_only_result_file must be specified.')
        usage()
        sys.exit()

    fptr = open(vina_score_only_result_file, 'r')
    lines = fptr.readlines()
    all_models = []
    num_models = 0
    for l in lines:
        if l.find("MODEL ")>-1:
            cur_model = {}               
            all_models.append(cur_model) #@@
            print("processing model %d:" %(all_models.index(cur_model)+1))
            num_models=len(all_models) 
        elif l.find("REMARK VINA RESULT:")>-1:  
            cur_model['energy'] = float(l.split()[3])   #-10.1
            print("    energy = ", cur_model['energy'])
        elif l.find('Affinity:')==0:
            cur_model['Affinity'] = float(l.split()[1]) #-10.29562
            print("    Affinity: ", cur_model['Affinity'])
        elif l.find("gauss 1")>-1:
            cur_model['gauss 1'] = float(l.split()[-1]) #119.70294
            print("    gauss 1: ", cur_model['gauss 1'])
        elif l.find("gauss 2")>-1:
            cur_model['gauss 2'] = float(l.split()[-1]) #2325.85768
            print("    gauss 2: ", cur_model['gauss 2'])
        elif l.find("repulsion")==4:
            cur_model['repulsion'] = float(l.split()[-1]) #2.42141
            print("    repulsion: ", cur_model['repulsion'])
        elif l.find("hydrophobic")==4:
            cur_model['hydrophobic'] = float(l.split()[-1])#78.61094
            print("    hydrophobic: ", cur_model['hydrophobic'])
        elif l.find("Hydrogen")==4:
            cur_model['Hydrogen'] = float(l.split()[-1]) #1.75464
            print("    Hydrogen: ", cur_model['Hydrogen'])
    num_models = len(all_models)
    for cur_model_index in range(num_models):
        cur_model = all_models[cur_model_index]
        if verbose: print("wt_dict[num_tors_div] = ", wt_dict['num_tors_div'])
        if verbose: print("vina_num_torsions = ", vina_num_torsions)
        num_tors_div = 1+ wt_dict['num_tors_div']* vina_num_torsions #@@ vina inored 2 of original 14 torsions
        if verbose: print("num_tors_div =", num_tors_div)
        print("Model %d calculated scores:" % ( cur_model_index + 1))
        #wt_dict['gauss 1'] = -0.035579
    	print("    weighted gauss 1 = %6.6f " % ( cur_model['gauss 1']*wt_dict['gauss 1'])) 
        #wt_dict['gauss 2'] = -0.005156
    	print("    weighted gauss 2 = %6.6f " % ( cur_model['gauss 2']*wt_dict['gauss 2']))
        #wt_dict['repulsion']= 0.840245
    	print("    weighted repulsion = %6.6f " % ( cur_model['repulsion']*wt_dict['repulsion']))
        #wt_dict['hydrophobic'] = -0.035069
    	print("    weighted hydrophobic = %6.6f " % ( cur_model['hydrophobic']*wt_dict['hydrophobic']))
        #wt_dict['Hydrogen'] = 0.0587439
    	print("    weighted Hydrogen = %6.6f " %( cur_model['Hydrogen']*wt_dict['Hydrogen']))
    	score = wt_dict['gauss 1']*cur_model['gauss 1'] + wt_dict['gauss 2']*cur_model['gauss 2'] + wt_dict['repulsion'] *cur_model ['repulsion']+ wt_dict['hydrophobic']*cur_model['hydrophobic'] + wt_dict['Hydrogen']*cur_model['Hydrogen'] 
    	print("    score without num_tors_div for model %d is %6.6f " % ( cur_model_index+1, score))   
        #num_tors_div = 1+ wt_dict['num_tors_div']* vina_num_torsions #@@ vina ignored 2 of original 14 torsions
        #print "    num_tors_div = % 6.6f " % ( num_tors_div )
        w = 0.1 * ( wt_dict['num_tors_div'] + 1)
        print("    w = ", w, " and vina_num_torsions =", vina_num_torsions)
        #print "    score/num_tors_div   is %6.6f " % ( score/num_tors_div)
        final_score = score/(1 + w * vina_num_torsions/5.0)
        print("    Final_score (score divided by (1+w *vina_num_torsions/5.0)) is %6.6f " % (final_score))
        # rh empirical: print "    alternative final score (score * 0.568137) is ", score * 0.568137 
    	cur_model['score'] = final_score
    	#print "score for model %d",  cur_model_index+1, " is ", cur_model['score']
        #print "model %d:" %(cur_model_index+1)
    	#print "   Final score is %6.6f " %( cur_model['score'])
        print("    Affinity  is %6.6f " %( cur_model['Affinity'])) 
    	#cur_model['score'] = final_score
    
    

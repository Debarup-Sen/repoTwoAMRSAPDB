#Sai
#OM Sai Ram

#Om Sai Ram
#Sai

import streamlit as lit
import subprocess as proc
from skbio.alignment import local_pairwise_align_protein as lalign, global_pairwise_align_protein as galign
from skbio import Protein
import pandas as pd
from io import StringIO

#lit.text(

lit.set_page_config(layout='wide', page_title="AMRSAP")

lit.write("""
# Welcome to the AMRSAP Sequence Alignment Toolbox!
*A toolbox for all alignment purposes.*

""")

dataset = lit.radio("Please select dataset against which to perform BLASTp (default: Non-clinical):",
                    ('Non-clinical Dataset', 'Clinical Dataset'))

tool = lit.radio(
    "Select an alignment option:",
    ('BLASTp (Basic Local Alignment Search against AMRSAPDB)', 'MUSCLE (Multiple Sequence Alignment)',
     'Needleman-Wunsch (Global Pairwise Alignment)', 'Smith-Waterman (Local Pairwise Alignment)'))

if 'BLASTp' in tool:
    query = lit.text_area('Enter your input protein sequence (in FASTA/multi-FASTA format/plain text sequence format/AMRSAPDB Acc. ID, e.g. AMRSAP1) here',
                          height=200).upper()
    query = query.replace(' ', '')
    file_query = lit.file_uploader("Or, you may upload file")#, label_visibility="collapsed")
    lit.markdown('<br>', unsafe_allow_html=True)
    outfmt = lit.radio(
        "Select an output format:",
        ('Default', 'Pairwise', 'Query-anchored showing identities', 
     'Query-anchored no identities', 'Flat query-anchored showing identities',
     'Flat query-anchored no identities', 'BLAST XML', 'Tabular', 'Tabular with comment lines',
     'Seqalign (Text ASN.1)', 'Seqalign (Binary ASN.1)', 'Comma-separated values',
     'BLAST archive (ASN.1)', 'Seqalign (JSON)')
        )
    outfmt = ('0' if 'Pairwise' in outfmt
              else '1' if 'Query-anchored showing identities' in outfmt
              else '2' if 'Query-anchored no identities' in outfmt
              else '3' if 'Flat query-anchored showing identities' in outfmt
              else '4' if 'Flat query-anchored no identities' in outfmt
              else '5' if 'BLAST XML' in outfmt
              else '6' if outfmt=='Tabular'
              else '7' if 'Tabular with comment lines' in outfmt
              else '8' if 'Seqalign (Text ASN.1)' in outfmt
              else '9' if 'Seqalign (Binary ASN.1)' in outfmt
              else '10' if 'Comma-separated values' in outfmt
              else '11' if 'BLAST archive (ASN.1)' in outfmt
              else '12' if 'Seqalign (JSON)' in outfmt
              else 'def' if 'Default' in outfmt
              else None)


    lit.text('Customization parameters & choices:')
    command = 'blastp -query blast_input.txt '
    task = 'blastp'
##    task = lit.radio("Type of task:",
##                     ('blastp',
##                      'blastp-fast (a faster version of blastp that uses a larger word size (6 vs 2-3))',
##                      'blastp-short (blastp optimized for queries shorter than 30 residues)'))
    if task: command += '-task '+task
    evalue = lit.text_input("Please enter e-value:")
    if evalue: command += ' -evalue '+evalue
    word_size = lit.text_input("Please enter word size:")
    if word_size: command += ' -word_size '+word_size
    gapopen = lit.text_input("Please enter gap opening penalty:")
    if gapopen:  command += ' -gapopen '+gapopen
    gapextend = lit.text_input("Please enter gap extension penalty:")
    if gapextend:  command += ' -gapextend '+gapextend
    matrix = lit.selectbox(
        "Please select the matrix:",
        ('',
         'BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90',
         'PAM30', 'PAM70', 'PAM250',
         'IDENTITY')
        )
    if matrix: command += ' -matrix '+matrix
    threshold = lit.text_input("Please enter minimum word score:")
    if threshold: command += '-threshold '+threshold
    num_alignments = lit.text_input("Please enter number of database sequences to show alignment for (Default 250):")
    if num_alignments: command += ' -num_alignments '+num_alignments
    qcov_hsp_perc = lit.text_input("Please enter percent query coverage per HSP:")
    if qcov_hsp_perc: command += ' -qcov_hsp_perc '+qcov_hsp_perc
    max_target_seqs = lit.text_input("Please enter maximum number of query to keep:")
    if max_target_seqs: command += ' -max_target_seqs '+max_target_seqs


    
###  Blast should reject "AMRSAP1" && add multi acc id 

    submit = lit.button('Submit')
  
    if (query or file_query) and submit:
        lit.info("Input has been successfully submitted. Please wait till processing is completed. Results will appear below.")
        
        if file_query:
            query = StringIO(file_query.getvalue().decode("utf-8")).read().upper()
            
        ######################################## Non-clinical dataset checks ########################################
        if dataset == "Non-clinical Dataset":
            if query and 'AMRSAP' in query and '>' not in query:
                check = [i for i in query.split('\n') if i!= '']
                for i in check:
                    if '_' in i and i.isalnum():
                        lit.info('If you are trying to enter AMRSAP Acc. ID, please make sure it is in correct format (e.g. AMRSAP1)')
                        query = None
                        break

            elif 'AMRSAP' not in query and '>' not in query and any(char.isdigit() for char in query):
                lit.info('Please re-check your input for invalid characters')
                query = None
            
                
            if query and 'AMRSAP' in query and '>' not in query:
                query = [i for i in query.replace(' ', '').split('\n') if i!='']
                with open('AMRSAPDB_Non_clinical_all_data.tsv') as f, open('blast_input.txt', 'w') as g:
                    for k in query:
                        f.seek(0,0)
                        l = ' '
                        while(True):
                            i = f.readline()
                            if i=='':
                                lit.error(f'The AMRSAPDB Acc. ID {k} does not match with our database. Please re-check')
                                query = None
                                break
                            j = i.split('\t')
                            if k in j[0]:
                                g.write('>'+k+'\n'+j[2]+'\n')                                
                                break
            elif query:
                open('blast_input.txt', 'w').write(query)

        ########################################## Clinical Dataset checks ##########################################
        if dataset == "Clinical Dataset":
            if query and 'AMRSAPFCS' in query and '>' not in query:
                check = [i for i in query.split('\n') if i!= '']
                for i in check:
                    if '_' in i and i.isalnum():
                        lit.info('If you are trying to enter AMRSAPFCS Acc. ID, please make sure it is in correct format (e.g. AMRSAPFCS1)')
                        query = None
                        break

            elif 'AMRSAPFCS' not in query and '>' not in query and any(char.isdigit() for char in query):
                lit.info('Please re-check your input for invalid characters')
                query = None
            
            if query and 'AMRSAPFCS' in query and '>' not in query:
                query = [i for i in query.replace(' ', '').split('\n') if i!='']
                with open('AMRSAPDB_Clinical_all_data.tsv') as f, open('blast_input.txt', 'w') as g:
                    for k in query:
                        f.seek(0,0)
                        l = ' '
                        while(True):
                            i = f.readline()
                            if i=='':
                                lit.error(f'The AMRSAPFCS Acc. ID {k} does not match with our database. Please re-check')
                                query = None
                                break
                            j = i.split('\t')
                            if k in j[1]:
                                g.write('>'+k+'\n'+j[3]+'\n')
                                break
            elif query:
                open('blast_input.txt', 'w').write(query)

        ############################################### Checks Done #################################################
        

        if query:
            command += (' -db AMRSAPDB' if dataset == "Non-clinical Dataset" else ' -db AMRSAPFCS')

            if outfmt == 'def':
                proc.Popen((command+' -out blast_output_def1 -outfmt 0').split())
                proc.run((command+' -out blast_output_def2 -outfmt 7').split())
            else:
                proc.run((command+' -out blast_output -outfmt '+outfmt).split())

   
            lit.info("Your output below: [Formats 7-13 show no output when no hits are found]")

            user_db = f'User selected databases: {dataset}\n'


            if outfmt == 'def':
                myFile = [i for i in open('blast_output_def2').readlines()]
                lit.text(''.join([i for i in myFile[:5] if '# Fields: ' not in i]))
                try:
                    headers = [i for i in myFile if '# Fields: ' in i][0].replace('# Fields: ','').split(',')
                    headers = (['Description', 'Source Organism', 'Target Organism'] +
                                headers[2:] +
                               ['AMRSAPDB Accession'])
                    data = [i.strip().split('\t') for i in myFile if '#' not in i]
                    my_file = (open('AMRSAPDB_Non_clinical_all_data.tsv') if dataset == "Non-clinical Dataset" else open('AMRSAPDB_Clinical_all_data.tsv')).readlines()
    ##                lit.write(data[i])
                    for i in range(len(data)):
                        for j in range(len(my_file)):
                            if str(data[i][1]) in my_file[j]:
                                line = my_file[j].split('\t')
                                accs = data[i][1]
                                data[i] = (
                                            [line[1]] +
                                            [line[4]] +
                                            [line[6].replace('"', '')] +
                                            data[i][2:] +
                                            [
                                                (
                                                    f'<a target="_blank" href="https://bblserver.org.in/tamrsar/amrsapdb/non-clinical/non-clinical-entry?id={accs}&db=non_clinical">{accs}</a>'
                                                    if dataset == "Non-clinical Dataset"
                                                    else f'<a target="_blank" href="https://bblserver.org.in/tamrsar/amrsapdb/clinical/clinical-entry?id={accs}&db=clinical">{accs}</a>'
                                                )
                                            ]
                                          )
                                break
                    myDF = pd.DataFrame(data, columns=headers)
                    del myFile, data, headers
                    lit.write("<div style='font-size: 14px; overflow-x: auto; overflow-y: auto'>" + myDF.to_html(escape=False, index=False) + "</div>", unsafe_allow_html=True)
                    lit.markdown('<br><br><u><b>*Alignments:*</b></u><br>', unsafe_allow_html=True)
                    myDF.to_csv('blast_output_def2', sep='\t')

                    output1 = ''.join(open('blast_output_def1').readlines()[18:])
                    crsr = 0
                    for i in range(len(output1)):
                        if '>' in output1[i]:
                            crsr = i
                            break
                    output1 = ''.join(output1[crsr:])
                    lit.text(output1)
                    open('blast_output_def1', 'w').write(output1)

                    open('blast_output', 'w').write(open('blast_output_def2').read().replace('user_database', user_db)
                                                    +'\n\nAlignments:\n'+
                                                    open('blast_output_def1').read())
                    lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.txt')
                except Exception as e:
                    lit.text("No hits found")
                    print(e)

            elif outfmt == '0':
                myFile = ''.join([i for i in open('blast_output').readlines()[18:]])
                lit.text(myFile.replace('user_database', user_db))
                myFile = open('blast_output').read().replace('user_database', user_db)
                open('blast_output', 'w').write(myFile)
                lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.txt')

            elif outfmt == '1':
                myFile = ''.join([i for i in open('blast_output').readlines()[18:]])
                lit.text(myFile.replace('user_database', user_db))
                myFile = open('blast_output').read().replace('user_database', user_db)
                open('blast_output', 'w').write(myFile)
                lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.txt')

            elif outfmt == '2':
                myFile = ''.join([i for i in open('blast_output').readlines()[18:]])
                lit.text(myFile.replace('user_database', user_db))
                myFile = open('blast_output').read().replace('user_database', user_db)
                open('blast_output', 'w').write(myFile)
                lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.txt')

            elif outfmt == '3':
                myFile = ''.join([i for i in open('blast_output').readlines()[18:]])
                lit.text(myFile.replace('user_database', user_db))
                myFile = open('blast_output').read().replace('user_database', user_db)
                open('blast_output', 'w').write(myFile)
                lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.txt')

            elif outfmt == '4':
                myFile = ''.join([i for i in open('blast_output').readlines()[18:]])
                lit.text(myFile.replace('user_database', user_db))
                myFile = open('blast_output').read().replace('user_database', user_db)
                open('blast_output', 'w').write(myFile)
                lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.txt')

            elif outfmt == '5':
                myFile = ''.join(open('blast_output').readlines())
                lit.text(myFile.replace('user_database', user_db))
                myFile = open('blast_output').read().replace('user_database', user_db)
                open('blast_output', 'w').write(myFile)
                lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.xml')
                
            elif outfmt == '6' or outfmt == '10':
                lit.text('(Please choose "Tabular with comment lines" to see column headers)')
                myFile = pd.DataFrame([(i.strip().split(',') if ',' in i else i.strip().split('\t')) for i in open('blast_output').readlines()])
                lit.table(myFile)
                if outfmt == '6':
                    myFile.to_csv('blast_output', sep='\t')
                    myFile = open('blast_output').read()
                    open('blast_output', 'w').write(user_db + myFile)
                    lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.tsv')
                elif outfmt == '10':
                    myFile.to_csv('blast_output')
                    myFile = open('blast_output').read()
                    open('blast_output', 'w').write(user_db + myFile)
                    lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.csv')

            elif outfmt == '7':
                myFile = [i for i in open('blast_output').readlines()]
                headers = [i for i in myFile if '# Fields: ' in i][0].replace('# Fields: ','').split(',')
                data = [i.strip().split('\t') for i in myFile if '#' not in i]
                myDF = pd.DataFrame(data, columns=headers)
                del myFile, data, headers
                lit.table(myDF)
                myDF.to_csv('blast_output', sep='\t')
                myFile = open('blast_output').read()
                open('blast_output', 'w').write(user_db + myFile)
                lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.csv')

            elif outfmt == '8':
                myFile = ''.join(open('blast_output').readlines())
                lit.text(myFile.replace('user_database', user_db))
                myFile = open('blast_output').read().replace('user_database', user_db)
                open('blast_output', 'w').write(myFile)
                lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.txt')

            elif outfmt == '9':
                lit.text('Binary output cannot be displayed in browser. Please download file to view output')
##                myFile = open('blast_output', 'rb').read()
##                open('blast_output', 'wb').write(user_db + myFile)
                lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out')

            elif outfmt == '11':
                myFile = ''.join(open('blast_output').readlines())
                lit.text(myFile.replace('user_database', user_db))
                myFile = open('blast_output').read().replace('user_database', user_db)
                open('blast_output', 'w').write(myFile)
                lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.asn')

            elif outfmt == '12':
                myFile = ''.join(open('blast_output').readlines())
                lit.text(myFile.replace('user_database', user_db))
                myFile = open('blast_output').read().replace('user_database', user_db)
                open('blast_output', 'w').write(myFile)
                lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out.json')
    ##        #########################################################
    ##        if  outfmt=='9':
    ##            lit.download_button("Download output file", open('blast_output', 'rb'), file_name='BLAST_out')
    ##        else:
    ##            lit.download_button("Download output file", open('blast_output'), file_name='BLAST_out')
    ##        #########################################################
    elif submit and not (query or file_query):
        lit.error("Please enter input sequence!")

if 'MUSCLE' in tool:
    multiseq = lit.text_area('Enter your input protein sequences (in FASTA format/multi-FASTA format/AMRSAPDB Acc. ID, e.g. AMRSAP111 one in each line) here:',
                            height=200).upper()
    file_query = lit.file_uploader("Or, you may upload file")#, label_visibility="collapsed")
    lit.markdown('<br>', unsafe_allow_html=True)
    
    command = 'muscle -in muscle_input.txt -out muscle_output '
    
    muscol1, muscol2 = lit.columns(2)
    with muscol1:
        maxiters = lit.text_input("Please enter maximum number of iterations:")
        if maxiters: command += ' -maxiters '+maxiters
        maxhours = lit.text_input("Please enter Maximum time to iterate in hours:")
        if maxhours: command += ' -maxhours '+maxhours
    with muscol2:
        outfmt = lit.radio(
                "Select an output format:",
                ('Default (FASTA)', 'HTML', 'GCG MSF', 
             'CLUSTALW', 'CLUSTALW with header'),
                horizontal = True)
        outfmt = ('-html ' if 'HTML' in outfmt
                  else '-msf ' if 'MSF' in outfmt
                  else '-clw ' if 'CLUSTALW' in outfmt
                  else '-clwstrict ' if 'header' in outfmt
                  else None)
        if outfmt:  command += outfmt

        smcol1, smcol2 = lit.columns(2)
        with smcol1:
            diags = lit.checkbox('Find diagonals?')
            if diags: command += ' -diags '
            stable = lit.checkbox('Output results in input order?')
            if stable: command += ' -stable '
        with smcol2:
            group = lit.checkbox('Group sequences by similarity?')
            if group: command += ' -group '


    submit = lit.button('Submit')
    if file_query:
        multiseq = StringIO(file_query.getvalue().decode("utf-8")).read().upper()
    if multiseq and submit:
        if len([i for i in multiseq.split('\n') if i!=''])>=1 and ('AMRSAP' or 'AMRSAPFCS') in multiseq:
            multiseq = [i for i in multiseq.replace(' ', '').split('\n') if i!='']
            with (open('AMRSAPDB_Non_clinical_all_data.tsv') if 'AMRSAP' in multiseq else open('AMRSAPDB_Clinical_all_data.tsv')) as f, open(r'muscle_input.txt', 'w') as g:
                for k in multiseq:
                    f.seek(0,0)
                    l = ' '
                    while(True):
                        i = f.readline()
                        if i=='':
                            lit.error(f'The AMRSAP Acc. ID {k} does not match with our database. Please re-check')
                            multiseq = None
                            break
                        j = i.split('\t')
                        if k in j[1]:
                            g.write('>'+k+'\n'+j[3]+'\n')
                            break
        else:
            open('muscle_input.txt', 'w').write(multiseq)
    
    if multiseq and submit:
        lit.info("Input has been successfully submitted. Please wait till processing is completed. Results will appear below.")
        proc.run((command+'').split())
        lit.info("Your output below (In case you do not see any output, please re-check your input for invalid characters or non-standard residues):")
        lit.text(open(r'muscle_output').read())
        lit.download_button("Download output file", open('muscle_output'), file_name='MUSCLE_out.txt')
    elif submit and not multiseq:
        lit.error("Please enter input sequence!")


if 'Needleman-Wunsch' in tool:
    lit.text("FASTA format, plain text sequence format supported.")
    query = myquery = lit.text_area('Enter your query protein sequence here:', height=200).upper()
    file_query = lit.file_uploader("Or, you may upload query file")#, label_visibility="collapsed")
    lit.markdown('<br>', unsafe_allow_html=True)
    subject = mysubject = lit.text_area('Enter your subject protein sequence here:', height=200).upper()
    file_subject = lit.file_uploader("Or, you may upload subject file")#, label_visibility="collapsed")
    gap_open_penalty = 0
    gap_extend_penalty = 0
    try:
        gap_open_penalty = float(lit.text_input("Please enter the gap open penalty: "))
    except:
        pass
    if not gap_open_penalty:
        gap_open_penalty = 11
    try:        
        gap_extend_penalty = float(lit.text_input("Please enter the gap extend penalty: "))
    except:
        pass
    if not gap_extend_penalty:
        gap_extend_penalty = 1
    submit = lit.button('Submit')
    if file_query:
        query = myquery = StringIO(file_query.getvalue().decode("utf-8")).read().upper()
    if file_subject:
        subject = mysubject = StringIO(file_subject.getvalue().decode("utf-8")).read().upper()
    if query and subject and submit:
        if '>' in query:
            query = ''.join(query.split('\n')[1:])
        if '>' in subject:
            subject = ''.join(subject.split('\n')[1:])
        if '\n' in query:
            query = query.replace('\n', '')
        if '\n' in subject:
            subject = subject.replace('\n', '')

        ######################################## Non-clinical dataset checks ########################################
        if dataset == "Non-clinical Dataset":
            if 'AMRSAP' not in query and query.isalpha() is False:
                lit.error("Some non-alphabet is present in the sequence. Please re-check!")
                query = None
            elif 'AMRSAP' in query and query.isalnum() is False:
                lit.error("Some unrecognized character is present in the Acc. ID. Please re-check!")
                query = None
            if 'AMRSAP' not in subject and subject.isalpha() is False:
                lit.error("Some non-alphabet is present in the sequence. Please re-check!")
                subject = None
            elif 'AMRSAP' in subject and subject.isalnum() is False:
                lit.error("Some unrecognized character is present in the Acc. ID. Please re-check!")
                subject = None
            if query and len([i for i in query.split('\n') if i!=''])==1 and 'AMRSAP' in query:
                with open('AMRSAPDB_Non_clinical_all_data.tsv') as f:
                    l = ' '
                    while(True):
                        i = f.readline()
                        if i=='':
                            lit.error(f'The AMRSAPDB Acc. ID {query} does not match with our database. Please re-check')
                            query = None
                            break
                        j = i.split('\t')
                        if query in j[1]:
                            query = j[3]
                            break
                        
            if subject and len([i for i in subject.split('\n') if i!=''])==1 and 'AMRSAP_' in subject:
                with open('AMRSAPDB_Non_clinical_all_data.tsv') as f:
                    l = ' '
                    while(True):
                        i = f.readline()
                        if i=='':
                            lit.error(f'The AMRSAPDB Acc. ID {subject} does not match with our database. Please re-check')
                            subject = None
                            break
                        j = i.split('\t')
                        if subject in j[1]:
                            subject = j[3]
                            break

        ######################################## Clinical dataset checks ########################################
        if dataset == "Clinical Dataset":
            if 'AMRSAPFCS' not in query and query.isalpha() is False:
                lit.error("Some non-alphabet is present in the sequence. Please re-check!")
                query = None
            elif 'AMRSAPFCS' in query and query.isalnum() is False:
                lit.error("Some unrecognized character is present in the Acc. ID. Please re-check!")
                query = None
            if 'AMRSAPFCS' not in subject and subject.isalpha() is False:
                lit.error("Some non-alphabet is present in the sequence. Please re-check!")
                subject = None
            elif 'AMRSAPFCS' in subject and subject.isalnum() is False:
                lit.error("Some unrecognized character is present in the Acc. ID. Please re-check!")
                subject = None
            if query and len([i for i in query.split('\n') if i!=''])==1 and 'AMRSAPFCS' in query:
                with open('AMRSAPDB_Clinical_all_data.tsv') as f:
                    l = ' '
                    while(True):
                        i = f.readline()
                        if i=='':
                            lit.error(f'The AMRSAPFCS Acc. ID {query} does not match with our database. Please re-check')
                            query = None
                            break
                        j = i.split('\t')
                        if query in j[1]:
                            query = j[3]
                            break
                        
            if subject and len([i for i in subject.split('\n') if i!=''])==1 and 'AMRSAPFCS' in subject:
                with open('AMRSAPDB_Clinical_all_data.tsv') as f:
                    l = ' '
                    while(True):
                        i = f.readline()
                        if i=='':
                            lit.error(f'The AMRSAPFCS Acc. ID {subject} does not match with our database. Please re-check')
                            subject = None
                            break
                        j = i.split('\t')
                        if subject in j[1]:
                            subject = j[3]
                            break

        ############################################### Checks Done #################################################

        if query and subject:
            lit.info("Input has been successfully submitted. Please wait till processing is completed. Results will appear below.")
            alignment, score, start_end_positions = galign(Protein(query.strip()), Protein(subject.strip()),
                                                           gap_open_penalty=gap_open_penalty, gap_extend_penalty=gap_extend_penalty)
            lit.info("Your output below:")
            lit.text(alignment)
            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.markdown('''<br>''', unsafe_allow_html=True)
            alignment.write(open('NWFile', 'w'))
            lines = [i for i in open('NWFile').readlines() if i!='']
            lit.text("Full alignment:")
            myquery = '>'+myquery+'\n'+lines[1]
            mysubject = '>'+mysubject+'\n'+lines[3]
            lit.text(myquery)
            lit.text(mysubject)
            lit.text("Score: "+str(score))
            my_write_file = 'AMRSAPDB Needleman-Wunsch Output:\n\nAlignment\n'+myquery+'\n'+mysubject+"\nScore: "+str(score)+'\n'
            open('NWFile', 'w').write(my_write_file)
            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.download_button("Download output file", open('NWFile'), file_name='NW_out')
    elif submit and (not query or not subject):
        lit.error("Please enter input sequence!")


if 'Smith-Waterman' in tool:
    lit.text("FASTA format, plain text sequence format supported.")
    query = myquery = lit.text_area('Enter your query protein sequence here', height=200).upper()
    file_query = lit.file_uploader("Or, you may upload query file")#, label_visibility="collapsed")
    lit.markdown('<br>', unsafe_allow_html=True)
    subject = mysubject = lit.text_area('Enter your subject protein sequence here', height=200).upper()
    file_subject = lit.file_uploader("Or, you may upload subject file")#, label_visibility="collapsed")
    gap_open_penalty = 0
    gap_extend_penalty = 0
    try:
        gap_open_penalty = float(lit.text_input("Please enter the gap open penalty: "))
    except:
        pass
    if not gap_open_penalty:
        gap_open_penalty = 11
    try:        
        gap_extend_penalty = float(lit.text_input("Please enter the gap extend penalty: "))
    except:
        pass
    if not gap_extend_penalty:
        gap_extend_penalty = 1
    submit = lit.button('Submit')
    if file_query:
        query = myquery = StringIO(file_query.getvalue().decode("utf-8")).read().upper()
    if file_subject:
        subject = mysubject = StringIO(file_subject.getvalue().decode("utf-8")).read().upper()
    if query and subject and submit:
        if '>' in query:
            query = ''.join(query.split('\n')[1:])
        if '>' in subject:
            subject = ''.join(subject.split('\n')[1:])
        if '\n' in query:
            query = query.replace('\n', '')
        if '\n' in subject:
            subject = subject.replace('\n', '')
        if 'AMRSAP' not in query and query.isalpha() is False:
            lit.error("Some non-alphabet is present in the sequence. Please re-check!")
            query = None
        elif 'AMRSAP' in query and query.replace('_', '').isalnum() is False:
            lit.error("Some unrecognized character is present in the Acc. ID. Please re-check!")
            query = None
        if 'AMRSAP' not in subject and subject.isalpha() is False:
            lit.error("Some non-alphabet is present in the sequence. Please re-check!")
            subject = None
        elif 'AMRSAP' in subject and subject.replace('_', '').isalnum() is False:
            lit.error("Some unrecognized character is present in the Acc. ID. Please re-check!")
            subject = None
        if query and len([i for i in query.split('\n') if i!=''])==1 and 'AMRSAP_' in query:
            with open('master_dataset.tsv') as f:
                l = ' '
                while(True):
                    i = f.readline()
                    if i=='':
                        lit.error(f'The AMRSAPDB Acc. ID {query} does not match with our database. Please re-check')
                        query = None
                        break
                    j = i.split('\t')
                    if query in j[1]:
                        query = j[3]
                        break
        if subject and len([i for i in subject.split('\n') if i!=''])==1 and 'AMRSAP_' in subject:
            with open('master_dataset.tsv') as f:
                l = ' '
                while(True):
                    i = f.readline()
                    if i=='':
                        lit.error(f'The AMRSAPDB Acc. ID {subject} does not match with our database. Please re-check')
                        subject = None
                        break
                    j = i.split('\t')
                    if subject in j[1]:
                        subject = j[3]
                        break
        if query and subject:
            lit.info("Input has been successfully submitted. Please wait till processing is completed. Results will appear below.")
            alignment, score, start_end_positions = lalign(Protein(query.strip()), Protein(subject.strip()),
                                                           gap_open_penalty=gap_open_penalty, gap_extend_penalty=gap_extend_penalty)
            lit.info("Your output below:")
            lit.text(alignment)
            alignment.write(open('SWFile', 'w'))
            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.markdown('''<br>''', unsafe_allow_html=True)
            alignment.write(open('SWFile', 'w'))
            lines = [i for i in open('SWFile').readlines() if i!='']
            lit.text("Full alignment:")
            myquery = '>'+myquery+'\n'+lines[1]
            mysubject = '>'+mysubject+'\n'+lines[3]
            lit.text(myquery)
            lit.text(mysubject)
            lit.text("Score: "+str(score))
            my_write_file = 'AMRSAPDB Smith-Waterman Output:\n\nAlignment\n'+myquery+'\n'+mysubject+"\nScore: "+str(score)+'\n'
            open('SWFile', 'w').write(my_write_file)
            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.markdown('''<br>''', unsafe_allow_html=True)
            lit.download_button("Download output file", open('SWFile'), file_name='SW_out')
    elif submit and not query:
        lit.error("Please enter input sequence!")

        

lit.write("*Thank you!*")

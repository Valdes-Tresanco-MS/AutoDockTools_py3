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
#  Modification date: 4/7/20 5:48                                                                  #
#                                                                                                  #
# ##################################################################################################

#############################################################################
#
# Author: Michel F. SANNER, Garrett MORRIS
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

PDBformat = {
    "HEADER": "HEADER    %.40s%.9s   %.4s",
    # continuation, repDate, idCode, rIdCode, rIdCode, rIdCode, rIdCode, rIdCode, rIdCode, rIdCode, rIdCode
    "OBSLTE": "OBSLTE%.2s  %.9s %.4s %.4s      %.4s %.4s %.4s %.4s %.4s %.4s %.4s ",
    # continuation, title
    "TITLE ": "TITLE   %.2s%60s",
    # continuation, idCode, comment
    "CAVEAT": "CAVEAT%.2s  %.4s %s    ",
    # continuation, compound
    "COMPND": "COMPND  %.2s%.60s",
    # continuation, srcName
    "SOURCE": "SOURCE  %.2s%.60s",
    # continuation, keywds
    "KEYWDS": "KEYWDS  %.2s%.60s",
    # continuation, technique
    "EXPDTA": "EXPDTA  %.2s%.60s",
    # modNum, continuation, modDate, modId, modType, record, record, record, record
    "REVDAT": "REVDAT %3d%.2s %.9s %.5s   %1d       %.6s %.6s %.6s %.6s",
    # continuation, sprsdeDate, idCode, sIdCode, sIdCode, sIdCode, sIdCode, sIdCode, sIdCode, sIdCode, sIdCode
    "SPRSDE": "SPRSDE%.2s  %.9s %.4s %.4s      %.4s %.4s %.4s %.4s %.4s %.4s %.4s ",
    # Journal
    "JRNL  ": "JRNL  %s      ",
    # remarkNum, empty
    "REMARK": "REMARK %3d %s",
    # idCode, chainID, seqBegin, insertBegin, seqEnd, insertEnd, database, dbAccession, dbIdCode, dbseqBegin, idbnsBeg,
    # dbseqEnd, dbinsEnd
    "DBREF ": "DBREF %.4s %.1s %4d %.1s%4d %.1s%s %s %s %5d %.1s%5d %.1s",
    # serNum, chainID, numRes, resName, resName, resName, resName, resName, resName, resName, resName, resName,
    # resName, resName, resName, resName
    "SEQRES": "SEQRES%2d  %.1s %4d %.3s  %.3s %.3s %.3s %.3s %.3s %.3s %.3s %.3s %.3s %.3s %.3s %.3s ",
    # idCode, resName, chainID, seqNum, iCode, stdRes, comment
    "MODRES": "MODRES%.4s %.3s %.1s %4d %.1s%.3s %s  ",
    # hetID, ChainID, seqNum, iCode, numHetAtoms, text
    "HET   ": "HET   %.3s %.1s  %4d%.1s%5d  %s     ",
    # continuation, hetID, text
    "HETNAM": "HETNAM%.2s  %.3s %s ",
    # continuation, hetID, hetSynonyms
    "HETSYN": "HETSYN%.2s  %.3s %.55s ",
    # compNum, hetID, continuation, asterisk, text
    "FORMUL": "FORMUL%2d  %.3s  %2d %.1s%s",
    # serNum, helixID, initResName, initChainID, initSeqNum, initICode, endResName, endChainID, endSeqNum, endICode,
    # helixClass, comment, length
    "HELIX ": "HELIX %3d %.3s %.3s %.1s %4d %.1s%.3s %.1s %4d %.1s%2d%s%5d ",
    # strand, sheetID, numStrands, initResName, initChainID, initSeqNum, initICode, endResName, endChainID, endSeqNum,
    # endICode, sense, curAtom, curResName, curChainId, curResSeq, curICode, prevAtom, prevResName, prevChainId,
    # prevResSeq, prevICode
    "SHEET ": "SHEET %3d %.3s %2d%.3s %.1s %4d%.1s%.3s %.1s %4d%.1s%2d%.4s %.3s%.1s %4d%.1s%.4s %.3s%.1s %4d%.1s",
    # seq, turnId, initResName, initChainId, initSeqNum, initICode, endResName, endChainId, endSeqNum, endICode,
    # comment
    "TURN  ": "TURN  %3d %.3s %.3s %.1s %4d%.1s%.3s %.1s %4d%.1s%s    ",
    # serNum, "CYS", chainID1, seqNum1, icode1, "CYS", chainID2, seqNum2, icode2, sym1, sym2
    "SSBOND": "SSBOND%3d %.3s %.1s %4d %.1s%.3s   %.1s %4d %.1s%.6s                       %.6s ",
    # name1, altLoc1, resName1, chainID1, resSeq1, iCode1, name2, altLoc2, resName2, chainID2, resSeq2, iCode2, sym1,
    # sym2
    "LINK  ": "LINK  %.4s      %.1s%.3s%.1s %4d%.1s%.4s               %.1s%.3s%.1s %4d%.1s%.6s  %.6s ",
    # name1, altLoc1, resName1, Chain1, resSeq1, ICode1, nameH, altLocH, ChainH, resSeqH, iCodeH, name2, altLoc2,
    # resName2, chainID2, resSeq2, iCode2, sym1, sym2
    "HYDBND": "HYDBND%.4s      %.1s%.3s%.1s %5d%.1s%.4s %.1s%.1s %5d%.1s%.4s %.1s%.3s%.1s %5d%.1s%.6s%.6s ",
    # atom1, altLoc1, resName1, chainID1, resSeq1, iCode1, atom2, altLoc2, resName2, chainID2, resSeq2, iCode2, sym1,
    # sym2
    "SLTBRG": "SLTBRG%.4s      %.1s%.3s%.1s %4d%.1s%.4s               %.1s%.3s%.1s %4d%.1s%.6s  %.6s ",
    # serNum, pep1, chainID1, seqNum1, icode1, pep2, chainID2, seqNum2, icode2, modNum, measure
    "CISPEP": "CISPEP%3d %.3s %.1s %4d %.1s%.3s   %.1s %4d %.1s%3d       %6.2f       ",
    # seqNum, siteID, numRes, resName1, chainID1, seq1, iCode1, resName2, chainID2, seq2, iCode2, resName3, chainID3,
    # seq3, iCode3, resName4, chainID4, seq4, iCode4
    "SITE  ": "SITE  %3d %.3s %2d %.3s %.1s %4d%.1s%.3s %.1s %4d%.1s%.3s %.1s %4d%.1s%.3s %.1s %4d%.1s",
    # a, b, c, alpha, beta, gamma, sGroup, z
    "CRYST1": "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s %4d",
    # o[n][1], o[n][2], o[n][3], t[n]
    "ORIGXn": "ORIGXn%10.6f    %10.6f%10.6f%10.5f     ",
    # o[n][1], o[n][2], o[n][3], t[n]
    "SCALEn": "SCALEn%10.6f    %10.6f%10.6f%10.5f     ",
    # serial, m[n][1], m[n][2], m[n][3], v[n], iGiven
    "MTRIXn": "MTRIXn%3d %10.6f%10.6f%10.6f%10.5f     %1d    ", "TVECT ": "TVECT %3d %10.5f%10.5f%10.5f%s",
    # serial
    "MODEL ": "MODEL %4d    ",
    # serial, name, altLoc, resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor, segID, element, charge
    "ATOM  ": "ATOM  %5d%.4s %.1s%.3s%.1s %4d%.1s%8.3f   %8.3f%8.3f%6.2f%6.2f%.4s      %.2s%.2s",
    # serial, name, altLoc, resName, chainID, resSeq, iCode, sigX, sigY, sigZ, sigOcc, sigTemp, segID, element, charge
    "SIGATM": "SIGATM%5d%.4s %.1s%.3s%.1s %4d%.1s%8.3f   %8.3f%8.3f%6.2f%6.2f%.4s      %.2s%.2s",
    # serial, name, altLoc, resName, chainID, resSeq, iCode, u[0][0], u[1][1], u[2][2], u[0][1], u[0][2], u[1][2],
    # segID, element, charge
    "ANISOU": "ANISOU%5d%.4s %.1s%.3s%.1s %4d%.1s%7d %7d%7d%7d%7d%7d%.4s  %.2s%.2s",
    # serial, name, altLoc, resName, chainID, resSeq, iCode, sig[1][1], sig[2][2], sig[3][3], sig[1][2], sig[1][3],
    # sig[2][3], segID, element, charge
    "SIGUIJ": "SIGUIJ%5d%.4s %.1s%.3s%.1s %4d%.1s%7d %7d%7d%7d%7d%7d%.4s  %.2s%.2s",
    # serial, resName, chainID, resSeq, iCode
    "TER   ": "TER   %5d%.3s      %.1s %4d%.1s",
    # serial, name, altLoc, resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor, segID, element, charge
    "HETATM": "HETATM%5d%.4s %.1s%.3s%.1s %4d%.1s%8.3f   %8.3f%8.3f%6.2f%6.2f%.4s      %.2s%.2s",
    "ENDMDL": "ENDMDL",
    # serial, serial, serial, serial, serial, serial, serial, serial, serial, serial, serial
    "CONECT": "CONECT%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d",
    # numRemark, 0, numHet, numHelix, numSheet, numTurn, numSite, numXform, numCoord, numTer, numConect, numSeq
    "MASTER": "MASTER%5d    %5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d"}


# PDBformat["JRNL"] = "JRNL  %.4s      %.2s%.51s " # "AUTH", continuation, authorList
# PDBformat["JRNL"] = "JRNL  %.4s      %.2s%s " # "TITL", continuation, title
# PDBformat["JRNL"] = "JRNL  %.4s      %.2s%.51s " # "EDIT", continuation, editorList
# PDBformat["JRNL"] = "JRNL  %.3s      %.15s   " # "REF", "TO BE PUBLISHED"
# PDBformat["JRNL"] = "JRNL  %.3s      %.2s%s %.2s  %s%s %4d " # "REF", continuation, pubName, "V.", volume, page, year
# PDBformat["JRNL"] = "JRNL  %.4s      %.2s%s " # "PUBL", continuation, pub
# PDBformat["JRNL"] = "JRNL  %.4s      %.4s                                                  " # "REFN", "0353"
# PDBformat["JRNL"] = "JRNL  %.4s      %.4s   %.6s %.2s  %.4s %s %.4s " # "REFN", "ASTM", astm, country, "ISBN", isbn,
# coden
# PDBformat["REMARK"] = "REMARK%.1s   %.9s %49d " # "1", "REFERENCE", refNum
# PDBformat["REMARK"] = "REMARK%.1s   %.4s  %.2s%.51s " # "1", "AUTH", continuation, authorList
# PDBformat["REMARK"] = "REMARK%.1s   %.4s  %.2s%s " # "1", "TITL", continuation, title
# PDBformat["REMARK"] = "REMARK%.1s   %.4s  %.2s%s " # "1", "EDIT", continuation, editorList
# PDBformat["REMARK"] = "REMARK%.1s   %.3s  %.15s   " # "1", "REF", "TO BE PUBLISHED"
# PDBformat["REMARK"] = "REMARK%.1s   %.3s  %.2s%s %.2s  %s%s %4d " # "1", "REF", continuation, pubName, "V.", volume,
# page, year
# PDBformat["REMARK"] = "REMARK%.1s   %.4s  %.2s%s " # "1", "PUBL", continuation, pub
# PDBformat["REMARK"] = "REMARK%.1s   %.4s  %.4s                                                  " # "1", "REFN",
# "0353"
# PDBformat["REMARK"] = "REMARK%.1s   %.4s  %.4s   %s %s  %.4s %s %.4s  " # "1", "REFN", "ASTM", astm, country, "ISBN",
# isbn, coden
# PDBformat["REMARK"] = "REMARK%.1s   %.11s %5.2f%.10s " # "2", "RESOLUTION.", resolution, "ANGSTROMS."
# PDBformat["REMARK"] = "REMARK%.1s   %.28s %s  " # "2", "RESOLUTION., comment
# PDBformat["REMARK"] = "REMARK%.1s   %.11s %s " # "2", "RESOLUTION.", comment
# FIXME
# PDBformat["SEQADV"] = "SEQADV%.4s %.3s %.1s %4d %.1s%s %s %.3s %5d %s " # idCode, resName, chainID, seqNum, iCode,
# database, dbIdCode, dbRes, dbSeq, conflict

PDBFormatConstr = {"SSBOND": [3, None, None, 4, None, None, None, 4, None, None, None],
                   "REVDAT": [3, None, None, None, 1, None, None, None, None], "MODEL ": [4],
                   "SCALEn": [10, 10, 10, 10], "COMPND": [None, None],
                   "SPRSDE": [None, None, None, None, None, None, None, None, None, None, None],
                   "CAVEAT": [None, None, None], "HETNAM": [None, None, None], "TVECT ": [3, 10, 10, 10, None],
                   "FORMUL": [2, None, 2, None, None], "SOURCE": [None, None],
                   "MODRES": [None, None, None, 4, None, None, None],
                   "HETATM": [5, None, None, None, None, 4, None, 8, 8, 8, 6, 6, None, None, None], "ENDMDL": [],
                   "ATOM  ": [5, None, None, None, None, 4, None, 8, 8, 8, 6, 6, None, None, None],
                   "CISPEP": [3, None, None, 4, None, None, None, 4, None, 3, 6], "CRYST1": [9, 9, 9, 7, 7, 7, None, 4],
                   "SLTBRG": [None, None, None, None, 4, None, None, None, None, None, 4, None, None, None],
                   "EXPDTA": [None, None],
                   "SIGUIJ": [5, None, None, None, None, 4, None, 7, 7, 7, 7, 7, 7, None, None, None],
                   "HET   ": [None, None, 4, None, 5, None],
                   "SHEET ": [3, None, 2, None, None, 4, None, None, None, 4, None, 2, None, None, None, 4, None, None,
                              None, None, 4, None],
                   "HYDBND": [None, None, None, None, 5, None, None, None, None, 5, None, None, None, None, None, 5,
                              None, None, None], "REMARK": [3, None], "TITLE ": [None, 60],
                   "ANISOU": [5, None, None, None, None, 4, None, 7, 7, 7, 7, 7, 7, None, None, None],
                   "SIGATM": [5, None, None, None, None, 4, None, 8, 8, 8, 6, 6, None, None, None],
                   "TURN  ": [3, None, None, None, 4, None, None, None, 4, None, None],
                   "SITE  ": [3, None, 2, None, None, 4, None, None, None, 4, None, None, None, 4, None, None, None, 4,
                              None], "HELIX ": [3, None, None, None, 4, None, None, None, 4, None, 2, None, 5],
                   "TER   ": [5, None, None, 4, None], "HEADER": [None, None, None],
                   "MASTER": [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5], "CONECT": [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
                   "HETSYN": [None, None, None],
                   "DBREF ": [None, None, 4, None, 4, None, None, None, None, 5, None, 5, None],
                   "MTRIXn": [3, 10, 10, 10, 10, 1], "ORIGXn": [10, 10, 10, 10], "JRNL  ": [None],
                   "SEQRES": [2, None, 4, None, None, None, None, None, None, None, None, None, None, None, None, None],
                   "LINK  ": [None, None, None, None, 4, None, None, None, None, None, 4, None, None, None],
                   "OBSLTE": [None, None, None, None, None, None, None, None, None, None, None], "KEYWDS": [None, None]}

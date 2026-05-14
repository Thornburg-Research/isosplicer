from pyLM import CME

import numpy as np
import pandas as pd
import json
import os

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import string

################################################################
def initSim(simVars):
    
    sim = CME.CMESimulation(name='splicing_toy')

    # Dummy reversible reaction to ensure the Gillespie algorithm timesteps remain smaller than the hybrid algorithm frequency
    species = ['DA','DB']
    sim.defineSpecies(species)
    sim.addParticles(species='DA', count=1000)
    sim.addParticles(species='DB', count=1000)
    sim.addReaction(reactant='DA', product='DB', rate=10)
    sim.addReaction(reactant='DB', product='DA', rate=10)

    # Initialize splicing species
    species = ['u1', 'u2', 'tri']
    sim.defineSpecies(species)
    
    sim.addParticles(species='u1', count=int(simVars['Values']['total_u1']*simVars['Values']['free_u1']))
    sim.addParticles(species='u2', count=int(simVars['Values']['total_u2']*simVars['Values']['free_u2']))
    sim.addParticles(species='tri', count=int(simVars['Values']['total_tri']*simVars['Values']['free_tri']))

    if simVars['Values']['pb'] != 0:
        for param in simVars.iterrows():
            if param[0].startswith('PTN_'):
                ptnID = param[0].split('PTN_')[1]
                print(ptnID)
                specPTN = [ptnID]
                sim.defineSpecies(specPTN)
                sim.addParticles(species=ptnID, count=int(simVars['Values'][param[0]]*simVars['Values'][f'free_{ptnID}']))
    
    return sim
################################################################
    

################################################################
def getGeneData(inDir, geneID, genes, mID=1):

    genes[geneID] = {}
    genes[geneID]['mID'] = mID

    
    geneFile = f'{inDir}{geneID}.xlsx'

    dataDF = pd.read_excel(geneFile)

    geneDF = dataDF.loc[dataDF['Feature'] == 'gene']
    geneStart = geneDF['Start'].values[0]
    print(geneStart)
    geneEnd = geneDF['End'].values[0]
    print(geneEnd)

    geneLen = geneEnd-geneStart

    genes[geneID]['length'] = int(geneLen)

    molec_mass = 337 #309 #g/mol/nucleotide
    density = 1.75*1000000 #g/m^3
    N_A = 6.023e23 #mol^-1
    
    R_H = ((3*molec_mass*int(geneLen))/(4*np.pi*6.023e23*density))**(1/3) # m
    mRNA_vol = ((4/3)*np.pi*(R_H**3))*1000 # L -> 1000 L/m^3
    genes[geneID]['mRNAV'] = mRNA_vol
    print('mRNA Volume: ', mRNA_vol)
    
    exonDF = dataDF.loc[dataDF['Feature'] == 'exon']

    if dataDF['Strand'].values[0] == '-':
        
        print('Reverse strand')

        exon5P = set(abs(exonDF['Start'].values - geneStart - geneLen))
        print(f'5 Prime Exons: {len(exon5P)}')
    
        # if len(exon5P)>26:
        #     raise Exception('Too many 5 prime splice sites for indexing')

        exon3P = set(abs(exonDF['End'].values - geneStart - geneLen))
        print(f'3 Prime Exons: {len(exon3P)}')
        
        # if len(exon3P)>26:
        #     raise Exception('Too many 3 prime splice sites for indexing')

    else:

        exon5P = set(exonDF['End'].values - geneStart)
        print(f'5 Prime Exons: {len(exon5P)}')
        
        # if len(exon5P)>26:
        #     raise Exception('Too many 5 prime splice sites for indexing')
    
        exon3P = set(exonDF['Start'].values - geneStart)
        print(f'3 Prime Exons: {len(exon3P)}')
    
        # if len(exon3P)>26:
        #     raise Exception('Too many 3 prime splice sites for indexing')

    
    allSites = list(sorted(exon5P|exon3P))

    print(len(allSites))
    print(allSites)

    idx5P = 0
    idx3P = 0
    exonTop = []
    
    for site in allSites:

        if site in exon3P:
            topID = index_to_letter(idx3P)
            print(site, '3P', topID)
            exonTop.append(topID)
            idx3P+=1
            # continue
        
        if site in exon5P:
            topID = index_to_letter(idx5P, uppercase=True)
            print(site, '5P', topID)
            exonTop.append(topID)
            idx5P+=1
            # continue


    startChecks = []
    endChecks = []
    
    for site in exonTop:
        if site[0].isupper():
            break
        startChecks.append(allSites[exonTop.index(site)])
    
    for site in reversed(exonTop):
        if site[0].islower():
            break
        # print(site)
        # print(exonTop.index(site),len(allSites))
        # endChecks.append(allSites[exonTop.index(site)])
        endChecks.append(allSites[-1])
    
    print(startChecks)
    print(endChecks)

    sites = []
    sitesIDX = {
        '3P':{},
        '5P':{}
    }

    for site in allSites:
        if (site not in startChecks) and (site not in endChecks):
            sites.append(site)
    
    print(len(sites))
    print(sites)


    idx5P = 0
    idx3P = 0
    allIdx = 0
    
    top = []
    
    for site in sites:

        if site in exon3P:
            topID = index_to_letter(idx3P)
            print(site, '3P', topID)
            top.append(topID)
            sitesIDX['3P'][site]=allIdx
            idx3P+=1
            allIdx+=1
            # continue
        
        if site in exon5P:
            topID = index_to_letter(idx5P, uppercase=True)
            print(site, '5P', topID)
            top.append(topID)
            sitesIDX['5P'][site]=allIdx
            idx5P+=1
            allIdx+=1
            # continue

    genes[geneID]['sites'] = sites
    genes[geneID]['siteidx'] = sitesIDX
    
    print(f'Topology: {top}')

    genes[geneID]['top'] = top

    return None
################################################################


################################################################
def getIsoforms(sim, genes, inDir, geneID, mID=1):

    geneFile = f'{inDir}{geneID}.xlsx'

    dataDF = pd.read_excel(geneFile)

    geneDF = dataDF.loc[dataDF['Feature'] == 'gene']
    geneStart = geneDF['Start'].values[0]
    print(geneStart)
    geneEnd = geneDF['End'].values[0]
    print(geneEnd)

    geneLen = geneEnd-geneStart

    genes[geneID]['length'] = int(geneLen)


    exonDF = dataDF.loc[dataDF['Feature'] == 'exon']

    if dataDF['Strand'].values[0] == '-':
        
        print('Reverse strand')

        exon5P = set(abs(exonDF['Start'].values - geneStart - geneLen))
        print(f'5 Prime Exons: {len(exon5P)}')
    
        # if len(exon5P)>26:
        #     raise Exception('Too many 5 prime splice sites for indexing')

        exon3P = set(abs(exonDF['End'].values - geneStart - geneLen))
        print(f'3 Prime Exons: {len(exon3P)}')
        
        # if len(exon3P)>26:
        #     raise Exception('Too many 3 prime splice sites for indexing')
    
    else:

        exon5P = set(exonDF['End'].values - geneStart)
        print(f'5 Prime Exons: {len(exon5P)}')
        
        # if len(exon5P)>26:
        #     raise Exception('Too many 5 prime splice sites for indexing')
    
        exon3P = set(exonDF['Start'].values - geneStart)
        print(f'3 Prime Exons: {len(exon3P)}')
    
        # if len(exon3P)>26:
        #     raise Exception('Too many 3 prime splice sites for indexing')

    
    print(exon3P)
    print(exon5P)
    allSites = list(sorted(exon5P|exon3P))

    print(len(allSites))
    print(allSites)

    idx5P = 0
    idx3P = 0
    exonTop = []
    
    for site in allSites:

        if site in exon3P:
            topID = index_to_letter(idx3P)
            print(site, '3P', topID)
            exonTop.append(topID)
            idx3P+=1
            # continue
        
        if site in exon5P:
            topID = index_to_letter(idx5P, uppercase=True)
            print(site, '5P', topID)
            exonTop.append(topID)
            idx5P+=1
            # continue


    startChecks = []
    endChecks = []
    
    for site in exonTop:
        if site[0].isupper():
            break
        startChecks.append(allSites[exonTop.index(site)])
    
    for site in reversed(exonTop):
        if site[0].islower():
            break
        # endChecks.append(allSites[exonTop.index(site)])
        endChecks.append(allSites[-1])
        
    
    print(startChecks)
    print(endChecks)
    

    mrnaIDs = set(exonDF['transcript_id'].values)

    mRNA = []

    for mrnaID in mrnaIDs:
        print(mrnaID)
        mrnaDF = exonDF.loc[exonDF['transcript_id']==mrnaID]
        # print(mrnaDF)
        if dataDF['Strand'].values[0] == '-':
            mrna5P = set(abs(mrnaDF['Start'].values - geneStart - geneLen))
            mrna3P = set(abs(mrnaDF['End'].values - geneStart - geneLen))
        else:
            mrna3P = set(abs(mrnaDF['Start'].values - geneStart))
            mrna5P = set(abs(mrnaDF['End'].values - geneStart))
        print(mrna5P)
        print(mrna3P)
    
        mrnaSites = list(sorted(mrna5P|mrna3P))
        print(mrnaSites)
    
        mrnaTop = []
        
        for site in mrnaSites:
    
            if (site not in startChecks) and (site not in endChecks):
    
                # siteIdx = genes[geneID]['sites'].index(site)
                if (len(mrnaTop) == 0) and (site not in genes[geneID]['siteidx']['5P']):
                    continue
                    
                siteIdx = getSpliceSiteIdx(genes, geneID, site, mrnaTop)

                topID = genes[geneID]['top'][siteIdx]
        
                mrnaTop.append(topID)
    
        print(mrnaTop)
        
        if len(mrnaTop) == 0:
            print(f'Single Exon mRNA: {mrnaID}')
            continue
            # for site in mrnaSites:
            #     siteIdx = getSpliceSiteIdx(genes, geneID, site, mrnaTop)
        
            #     topID = genes[geneID]['top'][siteIdx]
        
            #     mrnaTop = mrnaTop + topID
                
        else: 
            if mrnaTop[-1][0].isupper():
                mrnaTop = mrnaTop[:-1]
                print(mrnaTop)
        mRNA.append(mrnaTop)
    
    genes[geneID]['isoforms'] = mRNA

    for m in mRNA:
        mSpec = 'm_{}_'.format(geneID) +''.join(m)
        sim.defineSpecies([mSpec])

    getAllowedSplicePairs(genes, geneID)
        
    return None
################################################################


################################################################
def getSpliceSiteIdx(genes, geneID, site, mrnaTop):

    if len(mrnaTop) == 0:
        siteIdx = genes[geneID]['siteidx']['5P'][site]
    elif mrnaTop[-1][0].islower():
        siteIdx = genes[geneID]['siteidx']['5P'][site]
    else:
        siteIdx = genes[geneID]['siteidx']['3P'][site]

    return siteIdx
################################################################


################################################################
def getAllowedSplicePairs(genes, geneID):

    allowedPairs = []

    for mRNA in genes[geneID]['isoforms']:

        for i in range(len(mRNA)-1):

            spair = mRNA[i:i+2]

            if spair[0][0].isupper() and spair[1][0].islower():

                if spair not in allowedPairs:

                    allowedPairs.append(spair)

    genes[geneID]['allowedSplicePairs'] = allowedPairs

    return None
################################################################


################################################################
def trscRxns(sim, genes, simVars, coTrsc=False):

    for geneID, geneDict in genes.items():

        top = geneDict['top']
        mIDs = geneDict['mID']
    
        for mID in range(1,mIDs+1):
            trsc_act = 'TrscAct_{}_{}'.format(geneID,mID)
            sim.defineSpecies([trsc_act])
            
            if mID == 1:
                sim.addParticles(species=trsc_act, count=1)
                
        gID = 'g_{}'.format(geneID)
        
        species = [gID]
        sim.defineSpecies(species)
        
        sim.addParticles(species=gID, count=1)
        
        if not coTrsc:

            print('No Cotranscriptional Splicing')
        
            for mID in range(1,mIDs+1):
    
                unbound_premRNA = ['pm_{}_{}_{}'.format(geneID,mID,site) for site in top]
                sim.defineSpecies(unbound_premRNA)
    
                trsc_prods = unbound_premRNA
    
                trsc_act = 'TrscAct_{}_{}'.format(geneID,mID)
                if mID == mIDs:
                    next_mID = 1
                else:
                    next_mID = mID + 1
                trsc_next_act = 'TrscAct_{}_{}'.format(geneID,next_mID)
        #         trsc_prods.append(trsc_next_act)
    
                trsc_prods.append(gID)
    
                print(unbound_premRNA)
    
                sim.addReaction(reactant=tuple([gID,trsc_act]), product=tuple(trsc_prods), rate=simVars['Values']['trsc_init'])
                
        elif coTrsc:
            
            for mID in range(1,mIDs+1):
                
                unbound_premRNA = ['pm_{}_{}_{}'.format(geneID,mID,site) for site in top]
                sim.defineSpecies(unbound_premRNA)
                
                trsc_act = 'TrscAct_{}_{}'.format(geneID,mID)
                
                gID = 'g_{}'.format(geneID)
                
                trsc_prod = 'g_{}_{}'.format(geneID,mID)
                sim.defineSpecies([trsc_prod])
    
                sim.addReaction(reactant=tuple([gID,trsc_act]), product=trsc_prod, rate=simVars['Values']['trsc_init'])
            
    return None
################################################################
        

################################################################
def splicingRxnsUniform(sim, genes, simVars):

    for geneID, geneDict in genes.items():

        top = geneDict['top']
        mIDs = geneDict['mID']
    
        for mID in range(1,mIDs+1):
        
            u_spec = []
    
            for site in top:
    
                if site[0].isupper():
    
                    unbound = 'pm_{}_{}_{}'.format(geneID, mID,site)
                    print(unbound)
    
                    bound = 'pm_{}_{}_u{}'.format(geneID, mID,site)
                    print(bound)
                    sim.defineSpecies([bound])
    
                    reactants = [unbound, 'u1']
                    product = bound
    
                    sim.addReaction(reactant=tuple(reactants), product=product, rate=simVars['Values']['u1f']/simVars['Values']['AN']/simVars['Values']['V'])
                    sim.addReaction(reactant=product, product=tuple(reactants), rate=simVars['Values']['u1r'])
    
                    # tbound = 'pm_{}_{}_t{}'.format(geneID, mID,site)
                    # print(tbound)
                    # sim.defineSpecies([tbound])
    
                    # reactants = [bound, 'tri']
                    # product = tbound
    
                    # sim.addReaction(reactant=tuple(reactants), product=product, rate=simVars['triuf'])
                    # sim.addReaction(reactant=product, product=tuple(reactants), rate=simVars['triur'])
    
                    delete_site = 'D_{}_{}_{}'.format(geneID, mID,site)
                    sim.defineSpecies([delete_site])
    
                    sim.addReaction(reactant=tuple([delete_site,unbound]),product='',rate=simVars['Values']['instantaneous'])
                    sim.addReaction(reactant=tuple([delete_site,bound]),product='u1',rate=simVars['Values']['instantaneous'])
                    # sim.addReaction(reactant=tuple([delete_site,tbound]),product=tuple(['u1','tri']),rate=simVars['instantaneous'])
    
                if site[0].islower():
    
                    unbound = 'pm_{}_{}_{}'.format(geneID, mID,site)
                    print(unbound)
    
                #     bound = 'pm_{}_{}_u{}'.format(geneID, mID,site)
                #     print(bound)
                #     sim.defineSpecies([bound])
    
                #     reactants = [unbound, 'u2']
                #     product = bound
    
                #     sim.addReaction(reactant=tuple(reactants), product=product, rate=simVars['Values']['u2f']/simVars['Values']['AN']/simVars['Values']['V'])
                #     sim.addReaction(reactant=product, product=tuple(reactants), rate=simVars['Values']['u2r'])
    
                    delete_site = 'D_{}_{}_{}'.format(geneID, mID,site)
                    sim.defineSpecies([delete_site])
    
                    sim.addReaction(reactant=tuple([delete_site,unbound]),product='',rate=simVars['Values']['instantaneous'])
                #     sim.addReaction(reactant=tuple([delete_site,bound]),product='u2',rate=simVars['Values']['instantaneous'])
    
    
            for site1 in top:
    
                for site2 in top:
    
                    if site1[0].isupper() and site2[0].islower():
    
                        site1Idx = top.index(site1)
                        site2Idx = top.index(site2)
    
                        if site2Idx > site1Idx:

                            if 'allowedSplicePairs' in geneDict:

                                # pair = f"{site1}{site2}"
                                pair = [site1,site2]
    
                                if pair not in geneDict['allowedSplicePairs']:

                                    # print(f'Pair not allowed {pair}')
    
                                    continue
    
                            # tbound = 'pm_{}_{}_t{}'.format(geneID, mID,site1)
                            u1bound = 'pm_{}_{}_u{}'.format(geneID, mID, site1)
    
                            # u2bound = 'pm_{}_{}_u{}'.format(geneID, mID, site2)
                            site3p = 'pm_{}_{}_{}'.format(geneID, mID,site2)
    
                            Ecomplex = 'pm_{}_{}_E_{}{}'.format(geneID, mID, site1, site2)
    
                            if Ecomplex not in sim.species_id:
                                sim.defineSpecies([Ecomplex])

                            for sitePOS, siteIDX in geneDict['siteidx']['5P'].items():
                                if siteIDX == geneDict['top'].index(site1):
                                    position5P = int(sitePOS)

                            for sitePOS, siteIDX in geneDict['siteidx']['3P'].items():
                                if siteIDX == geneDict['top'].index(site2):
                                    position3P = int(sitePOS)

                            intronLength = int(position3P-position5P)
                            print(intronLength)

                            molec_mass = 337 #309 #g/mol/nucleotide
                            density = 1.75*1000000 #g/m^3
                            N_A = 6.023e23 #mol^-1
                            
                            R_H = ((3*molec_mass*int(intronLength))/(4*np.pi*6.023e23*density))**(1/3) # m
                            intronVolume = ((4/3)*np.pi*(R_H**3))*1000 # L -> 1000 L/m^3
                            # genes[geneID]['mRNAV'] = mRNA_vol
                            print('Intron Volume: ', intronVolume)
    
                            # sim.addReaction(reactant=tuple([u1bound,u2bound]), product=upair, rate=simVars['Values']['upsf']/simVars['Values']['AN']/geneDict['mRNAV'])
                            sim.addReaction(reactant=tuple([u1bound,site3p]), product=Ecomplex, rate=simVars['Values']['upsf']/simVars['Values']['AN']/intronVolume)
                            
                            sim.addReaction(reactant=Ecomplex, product=tuple([u1bound,site3p]), rate=simVars['Values']['upsr'])

                            Acomplex = 'pm_{}_{}_A_{}{}'.format(geneID, mID, site1, site2)

                            if Acomplex not in sim.species_id:
                                sim.defineSpecies([Acomplex])

                            sim.addReaction(reactant=tuple([Ecomplex, 'u2']), product=Acomplex, rate=simVars['Values']['u2f']/simVars['Values']['AN']/simVars['Values']['V'])

                            preBcomplex = 'pm_{}_{}_pB_{}{}'.format(geneID, mID, site1, site2)

                            if preBcomplex not in sim.species_id:
                                sim.defineSpecies([preBcomplex])

                            sim.addReaction(reactant=tuple([Acomplex,'tri']), product=preBcomplex, rate=simVars['Values']['triuf']/simVars['Values']['AN']/simVars['Values']['V'])
                            sim.addReaction(reactant=preBcomplex, product=tuple([Acomplex,'tri']), rate=simVars['Values']['triur'])
    
                            activatedBcomplex = 'pm_{}_{}_aB_{}{}'.format(geneID, mID, site1, site2)
    
                            products = [activatedBcomplex, 'u1']
    
                            if activatedBcomplex not in sim.species_id:
                                sim.defineSpecies([activatedBcomplex])
    
                            print(site1, site2)
                            if site2Idx-site1Idx>1:
                                print(top[site1Idx+1:site2Idx])
                                cutTop = top[site1Idx+1:site2Idx]
                                for csite in cutTop:
                                    products.append('D_{}_{}_{}'.format(geneID, mID,csite))
                            print()
    
                            sim.addReaction(reactant=preBcomplex, product=tuple(products), rate=simVars['Values']['act'])
    
                            spliced = 'pm_{}_{}_s_{}{}'.format(geneID, mID, site1, site2)
    
                            if spliced not in sim.species_id:
                                sim.defineSpecies([spliced])
    
                            sim.addReaction(reactant=activatedBcomplex, product=tuple([spliced,'u2','tri']), rate=simVars['Values']['splice'])

    return None
################################################################


################################################################
def splicingRxnsIndividualSites(sim, genes, simVars, paramsFile):

    for geneID, geneDict in genes.items():

        top = geneDict['top']
        mIDs = geneDict['mID']

        u1DF = pd.read_excel(paramsFile, sheet_name='U1', index_col=0)
        u2DF = pd.read_excel(paramsFile, sheet_name='U2', index_col=0)
        pairsDF = pd.read_excel(paramsFile, sheet_name='pairs', index_col=0)

        params_excel = pd.ExcelFile(paramsFile)
        addProtBinders = False
        if simVars['Values']['pb']>0:
            if 'prot' in params_excel.sheet_names:
                print('Protein binding parameters found')
                protDF = pd.read_excel(paramsFile, sheet_name='prot', index_col=0)
                addProtBinders = True
    
        for mID in range(1,mIDs+1):
        
            u_spec = []
    
            for site in top:
    
                if site[0].isupper():
    
                    unbound = 'pm_{}_{}_{}'.format(geneID, mID,site)
                    print(unbound)
    
                    bound = 'pm_{}_{}_u{}'.format(geneID, mID,site)
                    print(bound)
                    sim.defineSpecies([bound])
    
                    reactants = [unbound, 'u1']
                    product = bound
    
                    sim.addReaction(reactant=tuple(reactants), product=product, rate=u1DF[site]['u1f']/simVars['Values']['AN']/simVars['Values']['V'])
                    sim.addReaction(reactant=product, product=tuple(reactants), rate=u1DF[site]['u1r'])
    
                    # tbound = 'pm_{}_{}_t{}'.format(geneID, mID,site)
                    # print(tbound)
                    # sim.defineSpecies([tbound])
    
                    # reactants = [bound, 'tri']
                    # product = tbound
    
                    # sim.addReaction(reactant=tuple(reactants), product=product, rate=simVars['triuf'])
                    # sim.addReaction(reactant=product, product=tuple(reactants), rate=simVars['triur'])
    
                    delete_site = 'D_{}_{}_{}'.format(geneID, mID,site)
                    sim.defineSpecies([delete_site])
    
                    sim.addReaction(reactant=tuple([delete_site,unbound]),product='',rate=simVars['Values']['instantaneous'])
                    sim.addReaction(reactant=tuple([delete_site,bound]),product='u1',rate=simVars['Values']['instantaneous'])
                    # sim.addReaction(reactant=tuple([delete_site,tbound]),product=tuple(['u1','tri']),rate=simVars['instantaneous'])

                    if addProtBinders:
                        siteProtDF = protDF[site]
                        if siteProtDF.any():
                            for pbp in siteProtDF.keys():
                                if pbp.endswith('_on') and siteProtDF[pbp]!=0:
                                    ptnID = pbp.split('_on')[0]
                                    ptnBoundSite = 'pm_{}_{}_{}_{}'.format(geneID, mID,site,ptnID)
                                    sim.defineSpecies([ptnBoundSite])
                                    sim.addReaction(reactant=tuple([unbound, ptnID]), product=ptnBoundSite, rate=siteProtDF[ptnID+'_on']/simVars['Values']['AN']/simVars['Values']['V'])
                                    
                                    if siteProtDF[ptnID+'_off']>0:
                                        sim.addReaction(reactant=ptnBoundSite, product=tuple([unbound, ptnID]), rate=siteProtDF[ptnID+'_off'])
                                        
                                    if siteProtDF[ptnID+'_uf']>0:
                                        sim.addReaction(reactant=tuple([ptnBoundSite, 'u1']), product=tuple([bound, ptnID]), rate=siteProtDF[ptnID+'_uf']/simVars['Values']['AN']/simVars['Values']['V'])

                                        # if siteProtDF[ptnID+'_ur']>0:
                                        #     sim.addReaction(reactant=tuple([bound, ptnID]), product=tuple([ptnBoundSite, 'u1']), rate=siteProtDF[ptnID+'_ur'])
                                            
                
                if site[0].islower():
    
                    unbound = 'pm_{}_{}_{}'.format(geneID, mID,site)
                    print(unbound)
    
                    # bound = 'pm_{}_{}_u{}'.format(geneID, mID,site)
                    # print(bound)
                    # sim.defineSpecies([bound])
    
                    # reactants = [unbound, 'u2']
                    # product = bound
    
                    # sim.addReaction(reactant=tuple(reactants), product=product, rate=u2DF[site]['u2f']/simVars['Values']['AN']/simVars['Values']['V'])
                    # sim.addReaction(reactant=product, product=tuple(reactants), rate=u2DF[site]['u2r'])
    
                    delete_site = 'D_{}_{}_{}'.format(geneID, mID,site)
                    sim.defineSpecies([delete_site])
    
                    sim.addReaction(reactant=tuple([delete_site,unbound]),product='',rate=simVars['Values']['instantaneous'])
                    # sim.addReaction(reactant=tuple([delete_site,bound]),product='u2',rate=simVars['Values']['instantaneous'])

                    if addProtBinders:
                        siteProtDF = protDF[site]
                        if siteProtDF.any():
                            for pbp in siteProtDF.keys():
                                if pbp.endswith('_on') and siteProtDF[pbp]!=0:
                                    ptnID = pbp.split('_on')[0]
                                    ptnBoundSite = 'pm_{}_{}_{}_{}'.format(geneID, mID,site,ptnID)
                                    sim.defineSpecies([ptnBoundSite])
                                    sim.addReaction(reactant=tuple([unbound, ptnID]), product=ptnBoundSite, rate=siteProtDF[ptnID+'_on']/simVars['Values']['AN']/simVars['Values']['V'])
                                    
                                    if siteProtDF[ptnID+'_off']>0:
                                        sim.addReaction(reactant=ptnBoundSite, product=tuple([unbound, ptnID]), rate=siteProtDF[ptnID+'_off'])
                                        
                                    # if siteProtDF[ptnID+'_uf']>0:
                                    #     sim.addReaction(reactant=tuple([ptnBoundSite, 'u2']), product=tuple([bound, ptnID]), rate=siteProtDF[ptnID+'_uf']/simVars['Values']['AN']/simVars['Values']['V'])

                                        # if siteProtDF[ptnID+'_ur']>0:
                                        #     sim.addReaction(reactant=tuple([bound, ptnID]), product=tuple([ptnBoundSite, 'u2']), rate=siteProtDF[ptnID+'_ur'])
    
            for site1 in top:
    
                for site2 in top:
    
                    if site1[0].isupper() and site2[0].islower():
    
                        site1Idx = top.index(site1)
                        site2Idx = top.index(site2)
    
                        if site2Idx > site1Idx:

                            if 'allowedSplicePairs' in geneDict:

                                # pair = f"{site1}{site2}"
                                pair = [site1,site2]
    
                                if pair not in geneDict['allowedSplicePairs']:

                                    # print(f'Pair not allowed {pair}')
    
                                    continue

                            spair = f'{site1}{site2}'
    
                            # tbound = 'pm_{}_{}_t{}'.format(geneID, mID,site1)
                            u1bound = 'pm_{}_{}_u{}'.format(geneID, mID, site1)
    
                            # u2bound = 'pm_{}_{}_u{}'.format(geneID, mID, site2)
                            site3p = 'pm_{}_{}_{}'.format(geneID, mID,site2)
    
                            Ecomplex = 'pm_{}_{}_E_{}{}'.format(geneID, mID, site1, site2)
    
                            if Ecomplex not in sim.species_id:
                                sim.defineSpecies([Ecomplex])

                            for sitePOS, siteIDX in geneDict['siteidx']['5P'].items():
                                if siteIDX == geneDict['top'].index(site1):
                                    position5P = int(sitePOS)

                            for sitePOS, siteIDX in geneDict['siteidx']['3P'].items():
                                if siteIDX == geneDict['top'].index(site2):
                                    position3P = int(sitePOS)

                            intronLength = int(position3P-position5P)
                            print(intronLength)

                            molec_mass = 337 #309 #g/mol/nucleotide
                            density = 1.75*1000000 #g/m^3
                            N_A = 6.023e23 #mol^-1
                            
                            R_H = ((3*molec_mass*int(intronLength))/(4*np.pi*6.023e23*density))**(1/3) # m
                            intronVolume = ((4/3)*np.pi*(R_H**3))*1000 # L -> 1000 L/m^3
                            # genes[geneID]['mRNAV'] = mRNA_vol
                            print('Intron Volume: ', intronVolume)
    
                            # sim.addReaction(reactant=tuple([u1bound,u2bound]), product=upair, rate=simVars['Values']['upsf']/simVars['Values']['AN']/geneDict['mRNAV'])
                            # sim.addReaction(reactant=tuple([u1bound,u2bound]), product=upair, rate=simVars['Values']['upsf']/simVars['Values']['AN']/intronVolume)

                            sim.addReaction(reactant=tuple([u1bound,site3p]), product=Ecomplex, rate=pairsDF[spair]['upsf']/simVars['Values']['AN']/intronVolume)
                            sim.addReaction(reactant=Ecomplex, product=tuple([u1bound,site3p]), rate=pairsDF[spair]['upsr'])

                            Acomplex = 'pm_{}_{}_A_{}{}'.format(geneID, mID, site1, site2)

                            if Acomplex not in sim.species_id:
                                sim.defineSpecies([Acomplex])

                            sim.addReaction(reactant=tuple([Ecomplex, 'u2']), product=Acomplex, rate=u2DF[site2]['u2f']/simVars['Values']['AN']/simVars['Values']['V'])

                            if addProtBinders:
                                siteProtDF = protDF[site2]
                                if siteProtDF.any():
                                    # print(site2)
                                    for pbp in siteProtDF.keys():
                                        if pbp.endswith('_on') and siteProtDF[pbp]!=0:
                                            ptnID = pbp.split('_on')[0]
                                            ptnBoundSite = 'pm_{}_{}_{}_{}'.format(geneID, mID,site2,ptnID)

                                            EcomplexPTN = 'pm_{}_{}_E_{}{}_{}'.format(geneID, mID, site1, site2, ptnID)

                                            if EcomplexPTN not in sim.species_id:
                                                sim.defineSpecies([EcomplexPTN])

                                            sim.addReaction(reactant=tuple([Ecomplex, ptnID]), product=EcomplexPTN, rate=siteProtDF[ptnID+'_on']/simVars['Values']['AN']/simVars['Values']['V'])

                                            sim.addReaction(reactant=EcomplexPTN, product=tuple([Ecomplex, ptnID]), rate=siteProtDF[ptnID+'_off'])
                                            sim.addReaction(reactant=tuple([u1bound, ptnBoundSite]), product=EcomplexPTN, rate=pairsDF[spair]['upsf']/simVars['Values']['AN']/intronVolume)
                                            sim.addReaction(reactant=EcomplexPTN, product=tuple([u1bound, ptnBoundSite]), rate=pairsDF[spair]['upsr'])

                                            if siteProtDF[ptnID+'_uf']>0:
                                                sim.addReaction(reactant=tuple([EcomplexPTN, 'u2']), product=tuple([Acomplex, ptnID]), rate=siteProtDF[ptnID+'_uf']/simVars['Values']['AN']/simVars['Values']['V'])

                            #     ### TO DO: FIX CONDITIONAL IN CASE OF MULTIPLE PROTEIN BINDERS ###
                            #     else:
                            #         sim.addReaction(reactant=tuple([Ecomplex, 'u2']), product=Acomplex, rate=u2DF[site2]['u2f']/simVars['Values']['AN']/simVars['Values']['V'])


                            # else:
                            #     sim.addReaction(reactant=tuple([Ecomplex, 'u2']), product=Acomplex, rate=u2DF[site2]['u2f']/simVars['Values']['AN']/simVars['Values']['V'])

                            preBcomplex = 'pm_{}_{}_pB_{}{}'.format(geneID, mID, site1, site2)

                            if preBcomplex not in sim.species_id:
                                sim.defineSpecies([preBcomplex])

                            sim.addReaction(reactant=tuple([Acomplex,'tri']), product=preBcomplex, rate=pairsDF[spair]['triuf']/simVars['Values']['AN']/simVars['Values']['V'])
                            sim.addReaction(reactant=preBcomplex, product=tuple([Acomplex,'tri']), rate=pairsDF[spair]['triur'])
    
                            activatedBcomplex = 'pm_{}_{}_aB_{}{}'.format(geneID, mID, site1, site2)
    
                            products = [activatedBcomplex, 'u1']
    
                            if activatedBcomplex not in sim.species_id:
                                sim.defineSpecies([activatedBcomplex])
    
                            print(site1, site2)
                            if site2Idx-site1Idx>1:
                                print(top[site1Idx+1:site2Idx])
                                cutTop = top[site1Idx+1:site2Idx]
                                for csite in cutTop:
                                    products.append('D_{}_{}_{}'.format(geneID, mID,csite))
                            print()
    
                            sim.addReaction(reactant=preBcomplex, product=tuple(products), rate=pairsDF[spair]['act'])
    
                            spliced = 'pm_{}_{}_s_{}{}'.format(geneID, mID, site1, site2)
    
                            if spliced not in sim.species_id:
                                sim.defineSpecies([spliced])
    
                            sim.addReaction(reactant=activatedBcomplex, product=tuple([spliced,'u2','tri']), rate=pairsDF[spair]['splice'])

    return None
################################################################


################################################################
def exportSiteParameters(genes, simVars, outputDir='./', proteinBinders=0):

    for geneID, geneDict in genes.items():

        top = geneDict['top']

        u1Pf = {'Parameter':['u1f']}
        u1Pr = {'Parameter':['u1r']}
        u2Pf = {'Parameter':['u2f']}
        # u2Pr = {'Parameter':['u2r']}
        
        uSelPf = {'Parameter':['upsf']}
        uSelPr = {'Parameter':['upsr']}
        tPf = {'Parameter':['triuf']}
        tPr = {'Parameter':['triur']}
        actP = {'Parameter':['act']}
        splP = {'Parameter':['splice']}

        
        for site in top:

            if site[0].isupper():

                u1Pf[site] = [simVars['Values']['u1f']] #*simVars['Values']['AN']*simVars['Values']['V']]
                u1Pr[site] = [simVars['Values']['u1r']]

            if site[0].islower():

                u2Pf[site] = [simVars['Values']['u2f']] #*simVars['Values']['AN']*simVars['Values']['V']]
                # u2Pr[site] = [simVars['Values']['u2r']]


        for site1 in top:

            for site2 in top:

                if site1[0].isupper() and site2[0].islower():

                    site1Idx = top.index(site1)
                    site2Idx = top.index(site2)

                    if site2Idx > site1Idx:

                        if 'allowedSplicePairs' in geneDict:

                            # pair = f"{site1}{site2}"
                            pair = [site1,site2]

                            if pair not in geneDict['allowedSplicePairs']:

                                # print(f'Pair not allowed {pair}')

                                continue

                        spair = f'{site1}{site2}'

                        uSelPf[spair] = [simVars['Values']['upsf']] #*simVars['Values']['AN']*simVars['Values']['V']]
                        uSelPr[spair] = [simVars['Values']['upsr']]

                        tPf[spair] = [simVars['Values']['triuf']] #*simVars['Values']['AN']*simVars['Values']['V']]
                        tPr[spair] = [simVars['Values']['triur']]

                        actP[spair] = [simVars['Values']['act']]
                        
                        splP[spair] = [simVars['Values']['splice']]

        # print(u1Pf)
    
        u1fDF = pd.DataFrame(u1Pf)
        u1rDF = pd.DataFrame(u1Pr)
        u1DF = pd.concat([u1fDF,u1rDF], ignore_index=True)
        # u1DF.to_excel(f"{outputDir}combined.xlsx", sheet_name="U1", index=False)
    
        u2fDF = pd.DataFrame(u2Pf)
        # u2rDF = pd.DataFrame(u2Pr)
        # u2DF = pd.concat([u2fDF,u2rDF], ignore_index=True)
        u2DF = pd.concat([u2fDF], ignore_index=True)
        # u2DF.to_excel(f"{outputDir}combined.xlsx", sheet_name="U2", index=False)
    
        uSelFDF = pd.DataFrame(uSelPf)
        uSelRDF = pd.DataFrame(uSelPr)
        tfDF = pd.DataFrame(tPf)
        trDF = pd.DataFrame(tPr)
        actDF = pd.DataFrame(actP)
        splDF = pd.DataFrame(splP)
        pairDF = pd.concat([uSelFDF,uSelRDF,tfDF,trDF,actDF,splDF], ignore_index=True)
    
        if proteinBinders>0:
            # protF = {'Parameter':[f'protF{i}' for i in range(1,proteinBinders+1)]}
            # protR = {'Parameter':[f'protR{i}' for i in range(1,proteinBinders+1)]}
            simVars['Values']['pb'] = int(proteinBinders)

            protP = {'Parameter':[]}

            for i in range(proteinBinders):
                new_row = pd.DataFrame({'Values': [0]}, index=[f'PTN_prot{i+1}'])
                simVars = pd.concat([simVars, new_row])
                new_row = pd.DataFrame({'Values': [0.0]}, index=[f'free_prot{i+1}'])
                simVars = pd.concat([simVars, new_row])
                
                protP['Parameter'].append(f'prot{i+1}_on')
                protP['Parameter'].append(f'prot{i+1}_off')
                protP['Parameter'].append(f'prot{i+1}_uf')
                protP['Parameter'].append(f'prot{i+1}_blocks')
                
            for site in top:
                protP[site]=[0.0 for i in range(proteinBinders*4)]

            protDF = pd.DataFrame(protP)

        if not os.path.exists(f"{outputDir}{geneID}_parameters.xlsx"):
            exportGlobalParameters(genes, simVars, outputDir=outputDir)
    
        with pd.ExcelWriter(f"{outputDir}{geneID}_parameters.xlsx", engine='openpyxl', mode='a') as writer:
            u1DF.to_excel(writer, sheet_name="U1", index=False)
            u2DF.to_excel(writer, sheet_name="U2", index=False)
            pairDF.to_excel(writer, sheet_name="pairs", index=False)
            if proteinBinders>0:
                protDF.to_excel(writer, sheet_name='prot', index=False)

    return None
################################################################


################################################################
def mapIsoforms(sim, genes):
    
    isoforms = {}
    
    for geneID, geneDict in genes.items():
        
        mRNA = []
        determineMRNA(geneDict['top'], geneDict['top'], geneDict['sites'], mRNA)
        
        for m in mRNA:
            mSpec = 'm_{}_'.format(geneID) + m
            sim.defineSpecies([mSpec])
        
        geneDict['isoforms'] = mRNA
        
    return None
################################################################
        
    
################################################################
def determineMRNA(top, oTop, siteCoords, mRNA, sTop=''):
    
    remainingTops = []
    splicedTops = []
    mRNA_temp = []
    # print(sTop, top)
    
    for s1 in range(len(top)):
    
        for s2 in range(s1+1,len(top)):

            site1 = top[s1]
            site2 = top[s2]

            if site1.isupper() and site2.islower():

                site1Idx = top.find(site1)
                site2Idx = top.find(site2)

                if (site2Idx > site1Idx):

                    # p5Idx = oTop.index(site1)
                    # p3Idx = oTop.index(site2)

                    # genomicDistance = siteCoords[p3Idx] - siteCoords[p5Idx]

                    # if genomicDistance<=1200:

                    new_sTop = sTop+site1+site2
                    
                    splitTop = top.split(top[site1Idx:site2Idx+1])
    
                    remainingTop = splitTop[0] + splitTop[1]
    
                    remainingTops.append(remainingTop)
                    splicedTops.append(new_sTop)
                    
                    if len(remainingTop) == 0:
                        mRNA_temp.append(new_sTop)
                        
    if len(splicedTops) == 0:
        print('No New Splices: ', sTop, top)
        mRNA_temp.append(sTop)
                    
#     print(splicedTops)
#     print(remainingTops)
#     print()
    
    for m in mRNA_temp:
        
        # print('mRNA: ', m)
        removed_sites = []
        
        for i in range(len(m)-1):
            
            for j in range(len(m)-1):
                
                pair1 = m[i:i+2]
                pair2 = m[j:j+2]
#                 print(pair1, pair2)
                
                if pair1[0].isupper() and pair2[0].isupper():
                
                    p5_1 = oTop.index(pair1[0])
                    p3_1 = oTop.index(pair1[1])

                    p5_2 = oTop.index(pair2[0])
                    p3_2 = oTop.index(pair2[1])

                    if (p5_1>p5_2 and p3_1<p3_2):

                        removed_sites.append(pair1[0])
                        removed_sites.append(pair1[1])

                    elif (p5_2>p5_1 and p3_2<p3_1):

                        removed_sites.append(pair2[0])
                        removed_sites.append(pair2[1])
            
#         print(m)
#         print(removed_sites)
        for rSite in removed_sites:
            m = m.replace(rSite,"")
        new_m = "".join(sorted(m, key=lambda x: oTop.index(x)))
        # print(new_m)
        if new_m not in mRNA:
            mRNA.append(new_m)
#         print()
                  
    if len(remainingTops)>1:
        for i in range(len(remainingTops)):
            if len(remainingTops[i])>0:
                determineMRNA(remainingTops[i], oTop, siteCoords, mRNA, splicedTops[i])
                
################################################################


################################################################
def index_to_letter(index, uppercase=False):
    if index < 26:
        if uppercase:
            return chr(index + ord('A'))  # Convert to uppercase letter
        else:
            return chr(index + ord('a'))  # Convert to lowercase letter
    else:
        longIdx = int(index/26)
        if uppercase:
            return chr(index-26*(longIdx) + ord('A')) + str(longIdx+1)
        else:
            return chr(index-26*(longIdx) + ord('a')) + str(longIdx+1)
################################################################


################################################################
def setSimVars(varsFile=None, nucR=6e-6, cellR=9e-6):

    if varsFile is None:

        sf = 3
        # nucR = 6e-6
        # cellR = 9e-6
        cellV = 4/3*3.1415926535*((cellR)**3)*1000/1e-15
        nucV = int(4/3*3.1415926535*((nucR)**3)*1000/1e-15)
        u1 = 1661e-9*6.02e23*cellV*1e-15
        u1c = float(f"{u1:.{sf}g}")
        u2 = 1199e-9*6.02e23*cellV*1e-15
        u2c = float(f"{u2:.{sf}g}")
        u3 = 680e-9*6.02e23*cellV*1e-15
        u3c = float(f"{u3:.{sf}g}")

        simVars = {
            'u1f':[1.56e6], # Measurement
            'u1r':[0.015], # Measurement
            'u2f':[1e4], # round up order of magnitude: 0.15 min-1 (Lim, Hertel ATP dependent complex formation rate)/60 sec*min-1/(0.4 (40% Hela extract)*1200e-9 (M snrpb2 in Hela))
            # 'u2r':[0.062], # Measurement NOW ASSUMING IRREVERSIBLE
            'upsf':[1.0e3], # Order of magnitude calculated from measurement of 2 min intron half life (0.5 per min) assuming fully extended 70 nt intron as the "radius" for the reaction volume and 1 splice site pair for the concentration
            'upsr':[1.57], # Assuming a KD of 10^4; another option 1.57:Rate matching a measured dissociation rate of U1-70K
            'triuf':[1e4], # Inferred from measurement in yeast: 0.2 per min for a yeast nucleus 500 nm in diameter containing 100 trimers
            'triur':[1/(1*60)], # Measured - slow dwell time 1 min for u4
            'act':[1/(0.06*60)], # Assumed - corresponds to shortest dwell time measured for u4
            'splice':[0.067], # 15 second dwell time for U2+U5
            'trsc_init':[1.0], # Arbitrary until integrated into model including multiple genes, recommend values 0.01 - 100
            'trsc_elongation':[2000], # Elongation rate bp/min
            'total_u1':[u1c],
            'total_u2':[u2c],
            'total_tri':[u3c],
            'free_u1':[1.0],
            'free_u2':[1.0],
            'free_tri':[1.0],
            'V':[nucV],
            'pb':0
            }

        simVars = pd.DataFrame.from_dict(simVars).T
        simVars = simVars.rename(columns={0:'Values'})

    else:
        simVars = pd.read_excel(varsFile, sheet_name='global', index_col=0)
        # print(simVars)
    
    simVars['Values']['V'] = simVars['Values']['V']*1e-15 # L
    simVars.loc['AN'] = [6.02e23]
    
    simVars.loc['instantaneous'] = [1e10] # s-1

    return simVars
################################################################


################################################################
def exportGlobalParameters(genes, simVars, outputDir='./'):

    sv = simVars.copy()

    for geneID, geneDict in genes.items():
        
        sv['Values']['u1f'] = simVars['Values']['u1f'] #*simVars['Values']['AN']*simVars['Values']['V']
        sv['Values']['u2f'] = simVars['Values']['u2f'] #*simVars['Values']['AN']*simVars['Values']['V']
        sv['Values']['triuf'] = simVars['Values']['triuf'] #*simVars['Values']['AN']*simVars['Values']['V'] # N-1 s-1

        sv['Values']['V'] = simVars['Values']['V']/1e-15 # L

        sv.to_excel(f"{outputDir}{geneID}_parameters.xlsx",sheet_name='global')

    return None
################################################################


################################################################
def displayTopology(genes, geneID, exportFile=None, height=1.0, start=0, blockColor='red', siteColor='red', labelColor='black'):

    exonStart = 0
    intronStart = 0
    plot_data = []
    
    genePlotDat = {
        'block_widths': [],
        'gap_widths': [],
        'vline_positions': [],
        'vline_labels': []
    }
    plot_scaling_factor = 10/int(genes[geneID]['length'])
    
    for siteID in range(len(genes[geneID]['top'])):
        site = genes[geneID]['top'][siteID]
        # print(site)
        if site.isupper():
            
            if genes[geneID]['top'][siteID+1].islower():
                exonEnd = int(genes[geneID]['sites'][siteID])
                block_width = exonEnd - exonStart
                genePlotDat['block_widths'].append(int(block_width)*plot_scaling_factor)
                intronStart = int(exonEnd)
                
            elif genes[geneID]['top'][siteID+1].isupper():
                pass
    
            genePlotDat['vline_positions'].append(int(genes[geneID]['sites'][siteID])*plot_scaling_factor)
            genePlotDat['vline_labels'].append(site)
                
        elif site.islower():
            
            if site == genes[geneID]['top'][-1]:
                exonEnd = int(genes[geneID]['length'])
                block_width = exonEnd - int(genes[geneID]['sites'][siteID])
                genePlotDat['block_widths'].append(int(block_width)*plot_scaling_factor)
    
            if genes[geneID]['top'][siteID-1].isupper():
                intronEnd = int(genes[geneID]['sites'][siteID])
                gap_width = intronEnd - intronStart
                genePlotDat['gap_widths'].append(int(gap_width)*plot_scaling_factor)
                exonStart = int(intronEnd)
                
            elif genes[geneID]['top'][siteID-1].islower():
                pass
    
            genePlotDat['vline_positions'].append(int(genes[geneID]['sites'][siteID])*plot_scaling_factor)
            genePlotDat['vline_labels'].append(site)

    genePlotDat['block_labels'] = [i+1 for i in range(len(genePlotDat['block_widths']))]
    
    plot_data.append(genePlotDat)

    for isoform in genes[geneID]['isoforms']:
        exonStart = 0
        intronStart = 0
        current_position = 0
        mrnaPlotDat = {
            'block_widths': [],
            'vline_positions': [],
            'vline_labels': []
        }
        plot_scaling_factor = 10/int(genes[geneID]['length'])
    
        for isoID in range(len(isoform)):
            site = isoform[isoID]
            siteID = genes[geneID]['top'].index(site)
            # print(site)
            if site.isupper():
                
                if genes[geneID]['top'][siteID+1].islower():
                    exonEnd = int(genes[geneID]['sites'][siteID])
                    block_width = exonEnd - exonStart
                    mrnaPlotDat['block_widths'].append(int(block_width)*plot_scaling_factor)
                    intronStart = int(exonEnd)
                    current_position = current_position + int(block_width)
                    
                elif genes[geneID]['top'][siteID+1].isupper():
                    pass
        
                mrnaPlotDat['vline_positions'].append(int(current_position)*plot_scaling_factor)
                # print(int(genes[geneID]['sites'][siteID]-intronEnd))
                mrnaPlotDat['vline_labels'].append(site)
                    
            elif site.islower():
                
                if site == genes[geneID]['top'][-1]:
                    exonEnd = int(genes[geneID]['length'])
                    block_width = exonEnd - int(genes[geneID]['sites'][siteID])
                    mrnaPlotDat['block_widths'].append(int(block_width)*plot_scaling_factor)
        
                if genes[geneID]['top'][siteID-1].isupper():
                    intronEnd = int(genes[geneID]['sites'][siteID])
                    gap_width = intronEnd - intronStart
                    # geneDat['gap_widths'].append(int(gap_width)*plot_scaling_factor)
                    exonStart = int(intronEnd)
                    
                elif genes[geneID]['top'][siteID-1].islower():
                    pass
        
                mrnaPlotDat['vline_positions'].append(int(current_position)*plot_scaling_factor)
                mrnaPlotDat['vline_labels'].append(site)

        mrnaPlotDat['block_labels'] = [i+1 for i in range(len(mrnaPlotDat['block_widths']))]
    
        plot_data.append(mrnaPlotDat)

    # def draw_multiple_block_plots(plot_data, height=1, start=0):
    """
    Draws multiple vertically stacked block plots with labeled blocks and vertical lines.

    Args:
        plot_data (list of dict): Each dict must have:
            - block_widths (list of float)
            - gap_widths (list of float)
            - vline_positions (list of float, optional)
        height (float): Height of each rectangle.
        start (float): Starting x-position.
    """
    num_plots = len(plot_data)
    fig, axes = plt.subplots(nrows=num_plots, ncols=1, figsize=(12, 3 * num_plots), sharex=True)

    # If only one plot, axes is not a list
    if num_plots == 1:
        axes = [axes]

    for idx, data in enumerate(plot_data):
        ax = axes[idx]
        block_widths = data['block_widths']
        block_labels = data['block_labels']
        gap_widths = data.get('gap_widths', [])  # Optional now
        vline_positions = data.get('vline_positions', [])
        vline_labels = data['vline_labels']

        # If no gaps provided, use zeros
        if not gap_widths:
            gap_widths = [0] * (len(block_widths) - 1)

        if len(gap_widths) != len(block_widths) - 1:
            raise ValueError(f"Plot {idx}: Number of gaps must be one less than number of blocks.")

        total_length = sum(block_widths) + sum(gap_widths)

        # Total length of line
        total_length = sum(block_widths) + sum(gap_widths)

        # Draw base line
        ax.plot([start, start + total_length], [0, 0], color='black', linewidth=2, zorder=1)

        # Draw blocks with labels
        x = start
        for i, bw in enumerate(block_widths):
            rect = patches.Rectangle((x, -height/2), bw, height,
                                     linewidth=1, edgecolor='k', facecolor=blockColor, zorder=2)
            ax.add_patch(rect)

            # Centered number inside block
            # ax.text(x + bw / 2, 0, str(block_labels[i]), ha='center', va='center', fontsize=12, weight='bold', zorder=3)

            x += bw
            if i < len(gap_widths):
                x += gap_widths[i]

        # vline_height = 1
        # Draw vertical lines with alphabetic labels
        for i, xpos in enumerate(vline_positions):
            # label = string.ascii_uppercase[i % 26]
            label = vline_labels[i]
            if label.isupper():
                vls = 0.75
                offset=0.0
            elif label.islower():
                vls=-0.75
                offset=0.25
            ax.plot([xpos, xpos], [0, vls*height], color=siteColor, linestyle=':', linewidth=2)
            # ax.text(xpos+0.1, vls*(height / 2 + 0.3), label, color='red', ha='center', va='bottom', fontsize=10, weight='bold')
            ax.text(xpos, vls*height-offset, label, color=labelColor, ha='center', va='bottom', fontsize=16) # , weight='bold'
        
        ax.set_ylim(-height, height + 0.5)
        # ax.set_xlim(start - 1, start + total_length + 1)
        ax.set_aspect('equal')
        ax.axis('off')

    plt.tight_layout()
    plt.subplots_adjust(hspace=-0.8)
    
    if exportFile is not None:
        fig.set_dpi(300)
        plt.savefig(exportFile, transparent=True)
        
    plt.show()

    return None
################################################################

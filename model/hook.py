import lm
import numpy as np


################################################################
class SplicingSolver(lm.GillespieDSolver):

    ################################################################
    def __init__(self, pmap, genes, gparams = None, coTrsc = False):

        lm.GillespieDSolver.__init__(self)

        self.pmap = pmap

        self.genes = genes
        
        self.coTrsc = coTrsc

        self.gparams = gparams

        if gparams is None:
            self.trsc_elongation = 2000
        else:
            self.trsc_elongation = int(gparams['Values']['trsc_elongation'])

        # self.nextSpliceTime = 10.0

        if not self.coTrsc:
            print('Hook: No Cotranscriptional Splicing')

        return None
    ################################################################

    
    ################################################################
    def restart(self):
        
        for geneID, geneDict in self.genes.items():
            
            geneDict['trscStates'] = {}
            
            for mID in range(1,geneDict['mID']+1):
                
                geneDict['trscStates'][str(mID)] = {}
                geneDict['trscStates'][str(mID)]['position'] = 0
                geneDict['trscStates'][str(mID)]['lastSite'] = None
                geneDict['trscStates'][str(mID)]['lastTrscPos'] = None

        return None
    ################################################################

    
    ################################################################
    def hookSimulation(self, time):
        
        if time < 1.9:
            
            self.restart()

        # print('Time: ', time)

        self.counts = self.getSpeciesCountView()
#             print(time)

        if self.coTrsc:

            self.updateTranscription()
    
        self.resolveSplicing()

            # self.nextSpliceTime = self.nextSpliceTime + 1.0

        return 1
    ################################################################
    
    
    ################################################################
    def resolveSplicing(self):
        
        for geneID, geneDict in self.genes.items():
            
            mRNA = geneDict['isoforms']
            
            for mID in range(1,geneDict['mID']+1):

                for m in mRNA:

                    # print(m)

                    for i in range(len(m)-1):

                        spair = m[i:i+2]

                        if spair[0][0].isupper():

                            splice_state = 'pm_{}_{}_s_{}'.format(geneID, mID, ''.join(spair))
#                             print(geneID, mID, spair)
#                                 print(splice_state)

                            if self.counts[self.getIdx(splice_state)] == 0:
                                break

                        if i == int(len(m)-2):

                            mrnaID = 'm_{}_{}'.format(geneID, ''.join(m))

                            self.counts[self.getIdx(mrnaID)] = self.counts[self.getIdx(mrnaID)] + 1
                            
                            for m in mRNA:

                                for i in range(len(m)-1):

                                    spair = m[i:i+2]

                                    if spair[0][0].isupper():
                                        
                                        splice_state = 'pm_{}_{}_s_{}'.format(geneID, mID, ''.join(spair))
                                        
                                        self.counts[self.getIdx(splice_state)] = 0
                            
        return None
    ################################################################
    
    
    ################################################################
    def updateTranscription(self):
        
        for geneID, geneDict in self.genes.items():
            
            for mID in range(1,geneDict['mID']+1):
                
                trscID = 'g_{}_{}'.format(geneID,mID)
                
                if self.counts[self.getIdx(trscID)] > 0:

                    geneDict['trscStates'][str(mID)]['position'] = geneDict['trscStates'][str(mID)]['position'] + int(self.trsc_elongation/60)

                    # print(geneDict['trscStates'][str(mID)]['position'])
                        
                    currPos = 0

                    trscPos = []
            
                    for site in geneDict['sites']:
                        
                        if geneDict['trscStates'][str(mID)]['position'] > site:
                            
                            currPos = int(site)

                            if geneDict['trscStates'][str(mID)]['lastTrscPos'] is None:
                                trscPos.append(currPos)
                            else:
                                if currPos>geneDict['trscStates'][str(mID)]['lastTrscPos']:
                                    trscPos.append(currPos)

                    for trscSite in trscPos:
                            
                        if trscSite in geneDict['sites']:
    
                            spliceSiteIdxs = []
    
                            # print(currPos, trscSite)
    
                            if trscSite in geneDict['siteidx']['5P']:
                                spliceSiteIdxs.append(geneDict['siteidx']['5P'][trscSite])
                            if trscSite in geneDict['siteidx']['3P']:
                                spliceSiteIdxs.append(geneDict['siteidx']['3P'][trscSite])
    
                            for spliceSiteIdx in spliceSiteIdxs:
    
                                # spliceSiteIdx = geneDict['siteidx']['5P'].index(currPos)
        
                                spliceSite = geneDict['top'][spliceSiteIdx]
        
                                if spliceSite != geneDict['trscStates'][str(mID)]['lastSite']:
        
                                    siteID = 'pm_{}_{}_{}'.format(geneID,mID,spliceSite)
                                    
                                    # print('Creating Site: {}'.format(siteID))
        
                                    self.counts[self.getIdx(siteID)] = 1
        
                                    geneDict['trscStates'][str(mID)]['lastSite'] = spliceSite

                                    geneDict['trscStates'][str(mID)]['lastTrscPos'] = currPos
                
                                    if currPos == geneDict['sites'][0]:
                        
                                        if mID == geneDict['mID']:
                                
                                            nextmID = 1
                                    
                                        else:
                                            
                                            nextmID = int(mID + 1)
                                            
        #                                 print('Creating next transcription initiator')
                        
                                        nextTrscAct = 'TrscAct_{}_{}'.format(geneID,nextmID)
                            
        #                                 print(nextTrscAct)
                            
                                        # self.counts[self.getIdx(nextTrscAct)] = 1
        
                                        gPartID= 'g_{}'.format(geneID)
                                    
                                        self.counts[self.getIdx(gPartID)] = 1
                                        
                    
                    
                    if geneDict['trscStates'][str(mID)]['position'] > geneDict['length']:
                        
                        geneDict['trscStates'][str(mID)]['position'] = 0
                        geneDict['trscStates'][str(mID)]['lastSite'] = None
                        
                        self.counts[self.getIdx(trscID)] = 0

                        if self.gparams is not None:

                            for param in self.gparams.iterrows():
                                if param[0].startswith('free_'):
                                    particle = param[0].split('_')[1]

                                    print(f'Adding particles: {particle}')

                                    if f'total_{particle}' in self.gparams.index:

                                        self.counts[self.getIdx(particle)] = int(self.gparams['Values'][f'total_{particle}']) # self.counts[self.getIdx(particle)]

                                    if self.gparams['Values']['pb'] != 0:

                                        if f'PTN_{particle}' in self.gparams.index:
    
                                            self.counts[self.getIdx(particle)] = int(self.gparams['Values'][f'PTN_{particle}'])
        
        return None
    ################################################################
    
    
    ################################################################
    def getIdx(self, pID):

        return self.pmap.index(pID)
    ################################################################
    
################################################################

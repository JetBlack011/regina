# Written By Sam Ward in 2025 for the TXGL at UT Austin
# Based on 'block' construction from Angela Yuan and crossing edges from Andy Zhang
# 
#
# The purpose of this code is to build and a 3-Ball with a knot, specified by a PD code, 
# contained in the edges. It takes in a PD code and returns the edges within the 
# triangulation which contain the knot. We use the software package 'Regina', from 
# Ben Burton.

import regina

class Prism:
#this class initializes a prism with three tetrahedra and can return the various faces for gluing.

    def __init__(self, isR = False):
        self.isR = isR     #add orientation signifier

        # Create three tetrahedra
        self.t0 = triangulation.newTetrahedron()
        self.t1 = triangulation.newTetrahedron()
        self.t2 = triangulation.newTetrahedron()

        # Gluings for prism structure [0, 3, 1, 2]
        self.t0.join(2, self.t1, regina.Perm4(0, 3, 1, 2))  # Glue tet 0 to tet 1    
        self.t2.join(2, self.t1, regina.Perm4(0, 1, 2, 3))  # Glue tet 2 to tet 1

        self.tets = [self.t0.index(), self.t1.index(), self.t2.index()]
        #Assign w0,w1,w2 and v0,v1,v2 to vertices on the tets of the prism. They will be referenced from 0'th tetrahedra to 2nd and will contain the integers referencing the vertices (to make permutations easier).
        #w_0 - tet 0, vert 2 |
        self.w0 = [2]
        #w_1 - tet 0, vert 3 | tet 1, vert 2 |
        self.w1 = [3,2]
        #w_2 - tet 0, vert 1 | tet 1, vert 3 | tet 2, vert 3
        self.w2 = [1,3,3]
        #v_0 - tet 0, vert 0 | tet 1, vert 0 | tet 2, vert 0
        self.v0 = [0,0,0]
        #v_1 - tet 1, vert 1 | tet 2, vert 1 |
        self.v1 = [1,1]
        #v_2 - tet 2, vert 2 
        self.v2 = [2]
        
        if self.isR:
         #reassign vertices based on map from Prism to PrismR, so we end up swapping w and v, 0 for 1, and so-on
            self.w0 = [2] #self.v2
            self.w1 = [1,1] #self.v1
            self.w2 = [0,0,0] #self.v0
            self.v0 = [3,3,1] #self.w2 - swap places 0 and 2 for tetrahedron swap
            self.v1 = [2,3] #self.w1
            self.v2 = [2] #self.w0

            #reassign tet names
            self.t0, self.t2 = self.t2, self.t0


        #add top, bottom, diagonal, and vertical edge references, each with index in reference to the facet side - STORING THESE WILL CAUSE A SEGMENTATION ERROR use at own risk.
        #self.top_edges = [self.t0.edge(self.w1[0],self.w2[0]), self.t0.edge(self.w0[0],self.w2[0]), self.t0.edge(self.w0[0],self.w1[0])]
        #self.diag_edges = [self.t1.edge(self.v1[0],self.w2[1]), self.t2.edge(self.v0[2],self.w2[2]), self.t0.edge(self.v0[0],self.w1[0])]
        #self.bot_edges = [self.t2.edge(self.v1[1],self.v2[0]), self.t2.edge(self.v0[2],self.v2[0]), self.t2.edge(self.v0[2],self.v1[1])]

        #for vert edges, the index refers to the actual vertex, not a facet
        #self.vert_edges = [self.t0.edge(self.w0[0],self.v0[0]), self.t1.edge(self.w1[1],self.v1[0]), self.t2.edge(self.w2[2],self.v2[0])]

        #Edges by index
        self.top_edges = [self.t0.edge(self.w1[0],self.w2[0]).index(), self.t0.edge(self.w0[0],self.w2[0]).index(), self.t0.edge(self.w0[0],self.w1[0]).index()]
        self.diag_edges = [self.t1.edge(self.v1[0],self.w2[1]).index(), self.t2.edge(self.v0[2],self.w2[2]).index(), self.t0.edge(self.v0[0],self.w1[0]).index()]
        self.bot_edges = [self.t2.edge(self.v1[1],self.v2[0]).index(), self.t2.edge(self.v0[2],self.v2[0]).index(), self.t2.edge(self.v0[2],self.v1[1]).index()]

        #for vert edges, the index refers to the actual vertex, not a facet
        self.vert_edges = [self.t0.edge(self.w0[0],self.v0[0]).index(), self.t1.edge(self.w1[1],self.v1[0]).index(), self.t2.edge(self.w2[2],self.v2[0]).index()]


    def top_face(self):
        return self.t0.face(2, self.v0[0]) #[0, 1, 2]
    
    def bottom_face(self):
        return self.t2.face(2, self.w2[2]) #[0,1,2]
    
    def square_face(self, i):
    #this returns the triangles (as a list of integers, ordered from 0 to 2, w to v), then tetrahedra, then facets (as integers) making up the square faces. If True, then the side is blue, if False, then the side is red (this is for gluing purposes, it lets us know what way the diaginal is). The 'top' face and tet are 0 and 2. the 'bottom' is 1 and 3
        if i < 0 or i > 2:
            raise ValueError("INVALID INPUT?")
        if i == 2:
            #faces opposite w2
            s = [[self.w0[0], self.w1[0], self.v0[0]], [self.w1[1], self.v0[1], self.v1[0]], self.t0, self.t1, self.w2[0], self.w2[1], False]
            if self.isR: 
                s[-1] = True
            return s
        if i == 1:
            #faces opposite v1, w1
            s = [[self.w0[0], self.w2[0], self.v0[0]], [self.w2[-1], self.v0[-1], self.v2[-1]], self.t0, self.t2, self.w1[0], self.v1[1], True]
            if self.isR: 
                s[-1] = False
            return s
        else:#face opposites v0
            s = [[self.w1[1], self.w2[1], self.v1[0]], [self.w2[-1], self.v1[-1], self.v2[-1]], self.t1, self.t2, self.v0[1], self.v0[2], False]
            if self.isR: 
                s[-1] = True
            return s
    

    def glue_square(self, vert1, pris2, vert2):
        #this method uses the facets (of sides 0,1,2) to reference what sides are being glued
        #the inputs of prism must be a prism type, then using squareface we glue the triangles together.
        sqface1 = self.square_face(vert1)
        sqface2 = pris2.square_face(vert2)

        orient1 = sqface1[-1]
        orient2 = sqface2[-1]
        #if orient1 == orient2:
            #print('Error, cannot glue faces of the same orientation.')
        
        #else:
        tet1top = sqface1[2]  #get tetrahedra
        tet1bot = sqface1[3]
        tet2top = sqface2[2]
        tet2bot = sqface2[3]

            #Get the facets for the permutation
        facet1top = sqface1[4]
        facet1bot = sqface1[5]
        facet2top = sqface2[4]
        facet2bot = sqface2[5]

            #get permutations to glue first prism to second prism
        ptop = [0,0,0,0]
        ptop[facet1top] = facet2top
        ptop[sqface1[0][0]] = sqface2[0][0]
        ptop[sqface1[0][1]] = sqface2[0][1]
        ptop[sqface1[0][2]] = sqface2[0][2]
        permTop = regina.Perm4(ptop[0], ptop[1], ptop[2], ptop[3])

        pbot = [0,0,0,0]
        pbot[facet1bot] = facet2bot
        pbot[sqface1[1][0]] = sqface2[1][0]
        pbot[sqface1[1][1]] = sqface2[1][1]
        pbot[sqface1[1][2]] = sqface2[1][2]
        permBot = regina.Perm4(pbot[0], pbot[1], pbot[2], pbot[3])

        tet1top.join(facet1top, tet2top, permTop)  #glue square faces   
        tet1bot.join(facet1bot, tet2bot, permBot)



class Block:
    #this class creates a plus-shaped block and stores the indices of the prisms. It can return the prisms (i,j,k,l) for gluing, and the edges making up the embedded crossing. The first entry of i, j, k, l are the indexes of the block
    def __init__(self, index, strandi, strandj, strandk, strandl):
        self.index = index
        self.prisms = [Prism(True), Prism(False), Prism(True), Prism(False), Prism(True), Prism(False), Prism(True), Prism(False)]
        self.i = [strandi]
        self.j = [strandj]
        self.k = [strandk]
        self.l = [strandl]
    
        self.strands = [self.i, self.j, self.k , self.l]

    #make gluings - pin wheel!
    #from the diagram, if both arrows are pointing in, then the vertex is 2 and the side being glued is the one opposite the vertex 2 and so on.
    #gluings - (block index, facet):
        self.prisms[0].glue_square(0, self.prisms[1], 0)
        self.prisms[1].glue_square(1, self.prisms[2], 1)
        self.prisms[2].glue_square(0, self.prisms[3], 0)
        self.prisms[3].glue_square(1, self.prisms[4], 1)
        self.prisms[4].glue_square(0, self.prisms[5], 0)
        self.prisms[5].glue_square(1, self.prisms[6], 1)
        self.prisms[6].glue_square(0, self.prisms[7], 0)
        self.prisms[7].glue_square(1, self.prisms[0], 1)

    #assign tetrahedra to i,j,k,l
        self.i.extend([self.prisms[0], self.prisms[1]])
        self.l.extend([self.prisms[2], self.prisms[3]])
        self.k.extend([self.prisms[4], self.prisms[5]])
        self.j.extend([self.prisms[6], self.prisms[7]])

    #get edges - From Andy's crossing idea
    #top crossing, from l to j the top edges of (prism, facet): (7,0), (3,0)
    #bottom crossing, from i to k the top edges of (prism, facet): (0,0), (4,0)
        #these are the indices of the edges making up the crossing
        self.crossing = [self.prisms[0].diag_edges[0], self.prisms[4].diag_edges[0], self.prisms[7].top_edges[0], self.prisms[3].top_edges[0]]

    def glue_arm(self, strand, block2):
        #this function takes in two blocks and joins the arms of two blocks by strand
        if not any(strand in sl for sl in self.strands) or not any(strand in sl for sl in block2.strands):
           print('Error, strand not shared by one or both blocks.')
        else:
            self_ind = [[i, e.index(strand)] for i, e in enumerate(self.strands) if strand in e][0] #this is to retrieve the index of the list containing the faces of the blocks (by strand)
            block2_ind = [[i, e.index(strand)] for i, e in enumerate(block2.strands) if strand in e][0]
            #this gets a little complicated. it'll be the lower indexed prism to higher indexed prism - so 1 to 2
            self.strands[self_ind[0]][1].glue_square(2 ,block2.strands[block2_ind[0]][2], 2)
            self.strands[self_ind[0]][2].glue_square(2 ,block2.strands[block2_ind[0]][1], 2)
            print('Glued together strand: ', strand)



class Knot_builder:
    #this takes in a pd-code string of the form 'i j k l, i j k l, ...' and embeds the corresponding knot into a 3-ball
    def __init__(self, pdcode):
        self.blocks = []
        self.pdcode = pdcode
        self.knot_edges = [] #this stores the indices of edges making up the knot

        #parse pd_code input:
        pd_code_nums =  [list(map(int, block.split())) for block in self.pdcode.split(',')]
        num_blocks = len(pd_code_nums)
        for i in range(num_blocks):
            if len(pd_code_nums[i]) != 4:
                print('Error, pd code not in correct form.')
                break

        #make blocks and collect the block indices for each strand index, so strand 1 is shared by blocks 0 and 4, perhaps. 
        total_strands = (2 * num_blocks)
        strand_blocks = {n: [] for n in range(1, total_strands+1)} #populate with strands for references
        for n in range(num_blocks):
            i = pd_code_nums[n][0]
            j = pd_code_nums[n][1]
            k = pd_code_nums[n][2]
            l = pd_code_nums[n][3]
            block = Block(n,i,j,k,l)
            self.knot_edges.extend(block.crossing)
            self.blocks.append(block)

            strand_blocks[i].append(n) #assign block index to the strands it belongs to
            strand_blocks[j].append(n)
            strand_blocks[k].append(n)
            strand_blocks[l].append(n)
            
           
        #make gluings. Each strand appears twice, starting at 1
        for r in range(1, total_strands+1):
            n = strand_blocks[r][0]
            m = strand_blocks[r][1]
            block1 = self.blocks[n]
            block2 = self.blocks[m]
            block1.glue_arm(r, block2)

        #Fill it in/ cone it
        triangulation.finiteToIdeal()    


def main(): #example test
    global triangulation 
    triangulation = regina.Triangulation3()
    L = input("Please PD code in the form 'i j k l, i j k l, ...': ")
    K = Knot_builder(L)
    print(triangulation.isoSig())

main()

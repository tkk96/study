import matplotlib.pyplot as plt
import numpy as np
import os
from math import pi, sin, cos

def writeHeader(blockMeshDict):

    blockMeshDict.write("/*--------------------------------*- C++ -*----------------------------------*\ \n")
    blockMeshDict.write("| =========                 |                                                 | \n")
    blockMeshDict.write("| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n")
    blockMeshDict.write("|  \\    /   O peration     | Version:  4.x                                   | \n")
    blockMeshDict.write("|   \\  /    A nd           | Web:      www.OpenFOAM.org                      | \n")
    blockMeshDict.write("|    \\/     M anipulation  |                                                 | \n")
    blockMeshDict.write("\*---------------------------------------------------------------------------*/ \n")
    blockMeshDict.write("FoamFile \n")
    blockMeshDict.write("{ \n")
    blockMeshDict.write("    version     2.0; \n")
    blockMeshDict.write("    format      ascii; \n")
    blockMeshDict.write("    class       dictionary; \n")
    blockMeshDict.write("    object      blockMeshDict; \n")
    blockMeshDict.write("} \n")
    blockMeshDict.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n")
    blockMeshDict.write(" \n")
    blockMeshDict.write("convertToMeters 0.001; \n")
    blockMeshDict.write(" \n")

    return 0

def writeVertices(blockMeshDict, vertices):
    blockMeshDict.write("vertices \n")
    blockMeshDict.write("( \n")

    for i in range(len(vertices)):
        blockMeshDict.write("    ("+str(vertices[i][0])+" "+str(vertices[i][1])+" "+str(vertices[i][2])+")               //"+str(i)+" \n")

    blockMeshDict.write("); \n")
    blockMeshDict.write(" \n")

    return 0

def writeBlocks(blockMeshDict, blocks):
    blockMeshDict.write("blocks \n")
    blockMeshDict.write("( \n")

    for block in blocks:
        writeBlock(blockMeshDict, block)

    blockMeshDict.write("); \n\n")

    return 0

def writeBlock(blockMeshDict, block):
    blockMeshDict.write("\thex ("+str(block[0])+" "+str(block[1])+" "+str(block[2])+" "+str(block[3])+" "+str(block[0]+8)+" "+str(block[1]+8)+" "+str(block[2]+8)+" "+str(block[3]+8)+")\n")
    blockMeshDict.write("\t("+str(block[4])+" "+str(block[5])+" "+str(block[6])+")\n")
    blockMeshDict.write("\tsimpleGrading ("+str(block[7])+" 1 1)\n\n")

def writeEdges(blockMeshDict, edges):
    blockMeshDict.write("edges \n")
    blockMeshDict.write("( \n")

    for edge in edges:
        blockMeshDict.write("\tarc "+str(edge[0])+" "+str(edge[1])+" ("+str(edge[2])+" "+str(edge[3])+" "+str(edge[4])+")\n")

    blockMeshDict.write("); \n\n")

    return 0

def writePatches(blockMeshDict, patches, Lmap):
    blockMeshDict.write("boundary \n")
    blockMeshDict.write("( \n")

    blockMeshDict.write("    inlet\n")
    blockMeshDict.write("    {\n")
    blockMeshDict.write("        type mappedPatch;\n")
    blockMeshDict.write("        offset          (0 0 "+str(Lmap)+");\n")
    blockMeshDict.write("        sampleRegion    region0;\n")
    blockMeshDict.write("        sampleMode      nearestCell;\n")
    blockMeshDict.write("        samplePatch     none;\n")
    blockMeshDict.write("        faces\n")
    blockMeshDict.write("        (\n")
    for patch in patches[0]:
        blockMeshDict.write("\t\t\t("+str(patch[0])+" "+str(patch[1])+" "+str(patch[2])+" "+str(patch[3])+")\n")
    blockMeshDict.write("        );\n")
    blockMeshDict.write("    }\n")
    blockMeshDict.write("    outlet\n")
    blockMeshDict.write("    {\n")
    blockMeshDict.write("        type patch;\n")
    blockMeshDict.write("        faces\n")
    blockMeshDict.write("        (\n")
    for patch in patches[1]:
        blockMeshDict.write("\t\t\t("+str(patch[0])+" "+str(patch[1])+" "+str(patch[2])+" "+str(patch[3])+")\n")
    blockMeshDict.write("        );\n")
    blockMeshDict.write("    }\n")

    for i in range(len(centres)):
        blockMeshDict.write("    walls\n")
        blockMeshDict.write("    {\n")
        blockMeshDict.write("        type wall;\n")
        blockMeshDict.write("        faces\n")
        blockMeshDict.write("        (\n")
        for j in range(4):
            blockMeshDict.write("\t\t\t("+str(patches[2][j+i*4][0])+" "+str(patches[2][j+i*4][1])+" "+str(patches[2][j+i*4][2])+" "+str(patches[2][j+i*4][3])+")\n")
            
#            blockMeshDict.write("\t\t\t("+str(patches[2][i])+" "+str(patches[2][1])+" "+str(patches[2][2])+" "+str(patches[2][3])+")\n")
        blockMeshDict.write("        );\n")
        blockMeshDict.write("    }\n")




    blockMeshDict.write("); \n\n")

    return 0

def generateVertices(R_square, R_coolant, angles, Zb, Zt, centres):

    vertices = []

    for centre in centres:
        vertices.append([R_square + centre[0], -R_square + centre[1], Zb])
        vertices.append([-R_square + centre[0], -R_square + centre[1], Zb])
        vertices.append([-R_square + centre[0], R_square + centre[1], Zb])
        vertices.append([R_square + centre[0], R_square + centre[1], Zb])

        for angle in angles:
            vertices.append([R_coolant*cos(angle*pi/180.) + centre[0], R_coolant*sin(angle*pi/180.) + centre[1], Zb])

        for i in range(8):
            vertices.append([vertices[i][0]+ centre[0], vertices[i][1]+ centre[1], Zt])

    return vertices


def generateBlocks(Ns, Ni, Nz, grading, centres):

    blocks = []

    for i in range(len(centres)):
        blocks.append([ 1 + i*16, 0+ i*16, 3+ i*16, 2+ i*16, Ns, Ns, Nz, 1])
        blocks.append([ 0+ i*16, 4+ i*16, 7+ i*16, 3+ i*16, Ni, Ns, Nz, grading])
        blocks.append([ 3+ i*16, 7+ i*16, 6+ i*16, 2+ i*16, Ni, Ns, Nz, grading])
        blocks.append([ 2+ i*16, 6+ i*16, 5+ i*16, 1+ i*16, Ni, Ns, Nz, grading])
        blocks.append([ 1+ i*16, 5+ i*16, 4+ i*16, 0+ i*16, Ni, Ns, Nz, grading])

    return blocks

def generateEdges(curvature, eAngles, R_coolant, grading, centres):

    edges = []

    for i in range(len(centres)):
        edges.append([3+ i*16, 0+ i*16, curvature*cos(eAngles[0]*pi/180.)+ centres[i][0], curvature*sin(eAngles[0]*pi/180.)+ centres[i][1], Zb])
        edges.append([0+ i*16, 1+ i*16, curvature*cos(eAngles[1]*pi/180.)+ centres[i][0], curvature*sin(eAngles[1]*pi/180.)+ centres[i][1], Zb])
        edges.append([1+ i*16, 2+ i*16, curvature*cos(eAngles[2]*pi/180.)+ centres[i][0], curvature*sin(eAngles[2]*pi/180.)+ centres[i][1], Zb])
        edges.append([2+ i*16, 3+ i*16, curvature*cos(eAngles[3]*pi/180.)+ centres[i][0], curvature*sin(eAngles[3]*pi/180.)+ centres[i][1], Zb])

        edges.append([7+ i*16, 4+ i*16, R_coolant*cos(eAngles[0]*pi/180.)+ centres[i][0], R_coolant*sin(eAngles[0]*pi/180.)+ centres[i][1], Zb])
        edges.append([4+ i*16, 5+ i*16, R_coolant*cos(eAngles[1]*pi/180.)+ centres[i][0], R_coolant*sin(eAngles[1]*pi/180.)+ centres[i][1], Zb])
        edges.append([5+ i*16, 6+ i*16, R_coolant*cos(eAngles[2]*pi/180.)+ centres[i][0], R_coolant*sin(eAngles[2]*pi/180.)+ centres[i][1], Zb])
        edges.append([6+ i*16, 7+ i*16, R_coolant*cos(eAngles[3]*pi/180.)+ centres[i][0], R_coolant*sin(eAngles[3]*pi/180.)+ centres[i][1], Zb])  

        for j in range(8):
            edges.append([edges[j][0]+8+ i*16, edges[j][1]+8+ i*16, edges[j+ i*16][2], edges[j+ i*16][3], Zt])

    return edges

def generatePatches(blocks, edges, centres):

    patches = []

    inlet_patches=[]
    outlet_patches=[]
    wall_patches=[]

    for block in blocks:
        inlet_patches.append([block[0], block[1], block[2], block[3]])

        outlet_patches.append([block[0]+8, block[1]+8, block[2]+8, block[3]+8])

    for i in range(len(centres)):
        for j in range(4):
            wall_patches.append([edges[4+j+i*16][0], edges[4+j+i*16][1], edges[4+j+i*16][1]+8, edges[4+j+i*16][0]+8])


    patches.append(inlet_patches)
    patches.append(outlet_patches)
    patches.append(wall_patches)

    return patches
###############################################################################
#                      Input parameters
###############################################################################

D_pipe = 4.24    #pipe diameter

R_pipe = D_pipe/2.

Lmap = R_pipe*2.*4./1000.       #mapping distance

H_Mo = 8*D_pipe+760#4.24*8+760#214.313 + 214.312    #pipe length

grading = 0.1                      #grading

R_square = R_pipe*0.5

curvature = R_pipe * 0.67

Zb = 0.0

Zt = Zb + H_Mo

N=25                    #Mesh number in radial direction

Ns = int(N*0.4)

Ni = int(N*0.6)

Nz = 397*15               #Mesh number in axial direction

angles = [-45., -135., 135., 45.]

eAngles = [0., -90., 180., 90.]

centres =[]

centres.append([0., 0.])

'''
for i in range(6):
	centres.append([P*cos(60*i*pi/180),P*sin(60*i*pi/180)])
'''
###############################################################################
#                      generate parameters
###############################################################################

vertices = generateVertices(R_square, R_pipe, angles, Zb, Zt, centres)

blocks = generateBlocks(Ns, Ni, Nz, grading, centres)

edges = generateEdges(curvature, eAngles, R_pipe, grading, centres)

patches = generatePatches(blocks, edges, centres)

###############################################################################
#                      write parameters
###############################################################################

blockMeshDict = open('system/blockMeshDict','w')

writeHeader(blockMeshDict)

writeVertices(blockMeshDict, vertices)

writeBlocks(blockMeshDict, blocks)

writeEdges(blockMeshDict, edges)

writePatches(blockMeshDict, patches, Lmap)

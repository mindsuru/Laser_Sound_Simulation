import bpy
from scipy.spatial import Delaunay
import numpy as np
import bmesh
import json

def import_data_from_json(filename):
    with open(filename, 'r') as file:
        data = json.load(file)
    return np.array(data)

def runde_xy(punkt):
    x, y, z = punkt
    return round(x), round(y), z

def filter_points(points):
    # Check if the list is empty or has only one element
    if not points or len(points) == 1:
        return points

    # Initialize the filtered list with the first point
    filtered_points = [points[0]]

    # Iterate over the points starting from the second one
    for current_point in points[1:]:
        # Compare the x-value of the current point with the x-value of the last point in the filtered list
        if current_point[0] != filtered_points[-1][0]:
            filtered_points.append(current_point)

    return filtered_points

def create_grid_mesh(points, bm, radius):
    # Erstelle das Mesh
    
    vert_dict = {}
    for p in points:
        vert_dict[tuple(p[:2])] = bm.verts.new(p)

    # Erstelle Kanten basierend auf X- und Y-Nachbarn
    for x, y, z in points:
        p = (x, y)
        for dx, dy in [(1, 0), (0, 1)]:
            neighbor = (x + dx, y + dy)
            if neighbor in vert_dict:
                try:
                    bm.edges.new([vert_dict[p], vert_dict[neighbor]])
                except ValueError:
                    # Kante existiert bereits
                    continue

    # Neue Logik zum Finden der äußersten Punkte
    max_x_pro_y = {}
    min_x_pro_y = {}
    outerPoints = []
    for punkt in points:
        x, y, z = punkt
        if y not in max_x_pro_y or x > max_x_pro_y[y][0]:
            max_x_pro_y[y] = punkt
        if y not in min_x_pro_y or x < min_x_pro_y[y][0]:
            min_x_pro_y[y] = punkt

    sortierte_y_werte = sorted(max_x_pro_y.keys(), reverse=True)

    aeußerste_punkte = [max_x_pro_y[y] for y in sortierte_y_werte]

    for y in reversed(sortierte_y_werte):
        if not np.array_equal(min_x_pro_y[y], max_x_pro_y[y]):
            aeußerste_punkte.append(min_x_pro_y[y])

    # Verbinde die äußersten Punkte mit Kanten
    for i in range(len(aeußerste_punkte)):
        p1 = tuple(aeußerste_punkte[i][:2])
        p2 = tuple(aeußerste_punkte[(i + 1) % len(aeußerste_punkte)][:2])
        if p1 in vert_dict and p2 in vert_dict:
            try:
                bm.edges.new([vert_dict[p1], vert_dict[p2]])
            except ValueError:
                outerPoints.append(p1)
                continue
            
    filteredOuterPoints = []
    filteredOuterPoints = filter_points(outerPoints)
    facesUL = []
    facesUR = []
    # Hinzufügen der neuen Kanten für diagonale Verbindungen innerhalb von Quadraten
    for x, y, z in points:
        if x < 0 and y > 0 or x > 0 and y < 0: 
            for dx, dy in [(-1, -1), (1, 1)]:
                neighbor = (x + dx, y + dy)
                if neighbor in vert_dict:
                    try:
                        bm.edges.new([vert_dict[(x, y)], vert_dict[neighbor]])
                        facesUL.append(((x, y), neighbor))
                    except ValueError:
                        
                        continue
        if (x < 0 and y < 0 or x > 0 and y > 0):
            for dx, dy in [(1, -1), (-1, 1)]:
                neighbor = (x + dx, y + dy)
                if neighbor in vert_dict:
                    try:
                        bm.edges.new([vert_dict[(x, y)], vert_dict[neighbor]])
                        facesUR.append(((x, y), neighbor))
                    except ValueError:
                        # Kante existiert bereits
                        continue
    
    leftDown = []
    leftUp = []
    rightUp = []
    rightDown = []
    filteredOuterPoints.append((-radius, 0))
    for p in filteredOuterPoints:
        if p[0] < 0 and p[1] <= 0:
            leftDown.append(p)
            leftUp.append((p[0], -p[1]))
            rightUp.append((-p[0], -p[1]))
            rightDown.append((-p[0], p[1]))

    peakPoint = rightUp[len(rightUp) - 1]
    peakCorner = rightUp[len(rightUp) - 2]
    
    for i in range(peakPoint[1] + 1, int(peakCorner[1])):
        try:
            print(peakPoint[0], peakPoint[1], peakCorner[0], i)
            bm.edges.new([vert_dict[peakPoint[1], peakPoint[0]], vert_dict[i, peakCorner[0]]])
            bm.edges.new([vert_dict[-peakPoint[1], peakPoint[0]], vert_dict[-i, peakCorner[0]]])
            bm.edges.new([vert_dict[peakPoint[1], -peakPoint[0]], vert_dict[i, -peakCorner[0]]])
            bm.edges.new([vert_dict[-peakPoint[1], -peakPoint[0]], vert_dict[-i, -peakCorner[0]]])
            bm.edges.new([vert_dict[-peakPoint[0], peakPoint[1]], vert_dict[-peakCorner[0], -i]])
            bm.edges.new([vert_dict[-peakPoint[0], peakPoint[1]], vert_dict[-peakCorner[0], i]])
            bm.edges.new([vert_dict[peakPoint[0], peakPoint[1]], vert_dict[peakCorner[0], i]])
            bm.edges.new([vert_dict[peakPoint[0], peakPoint[1]], vert_dict[peakCorner[0], -i]])
        except:
            continue
        
    allP = []
    allP.append(leftDown)
    allP.append(leftUp)
    allP.append(rightUp)
    allP.append(rightDown)

    for lists in allP: 
        for i in range(len(lists) - 1):
            try:
                # Erstellen einer neuen Kante zwischen den Punkten
                bm.edges.new([vert_dict[lists[i]], vert_dict[lists[i + 1]]])
            except Exception as e:
                # Ausgabe der Fehlermeldung
                print('error:', e)
                

    try:
        bm.edges.new([vert_dict[(-1, 0)], vert_dict[(0, 1)]])
        bm.edges.new([vert_dict[(0, 1)], vert_dict[(1, 0)]])
        bm.edges.new([vert_dict[(1, 0)], vert_dict[(0, -1)]])
        bm.edges.new([vert_dict[(0, -1)], vert_dict[(-1, 0)]])
    except:
        print('error')

    bm = add_diagonal_facesUL(bm, vert_dict, facesUL)
    bm = add_diagonal_facesUR(bm, vert_dict, facesUR)
    bm = add_fucking_facesCenter(bm, vert_dict)
    #bm.faces.new([vert_dict[(0, 0)], vert_dict[(0, 1)], vert_dict[(1, 0)]])
    return bm

def add_fucking_facesCenter(bm, vert_dict):
    
    bm.faces.new([vert_dict[(0, 0)], vert_dict[(-1, 0)], vert_dict[(0, 1)]])
    bm.faces.new([vert_dict[(-1, 0)], vert_dict[(-1, 1)], vert_dict[(0, 1)]])
    bm.faces.new([vert_dict[(0, 0)], vert_dict[(0, 1)], vert_dict[(1, 0)]])
    bm.faces.new([vert_dict[(0, 1)], vert_dict[(1, 1)], vert_dict[(1, 0)]])
    bm.faces.new([vert_dict[(0, 0)], vert_dict[(1, 0)], vert_dict[(0, -1)]])
    bm.faces.new([vert_dict[(1, 0)], vert_dict[(1, -1)], vert_dict[(0, -1)]])
    bm.faces.new([vert_dict[(0, 0)], vert_dict[(0, -1)], vert_dict[(-1, 0)]])
    bm.faces.new([vert_dict[(-1, -1)], vert_dict[(-1, 0)], vert_dict[(0, -1)]])
    
    return bm
    
    
def add_diagonal_facesUR(bm, vert_dict, facesUR):
    for (x, y), (nx, ny) in facesUR:
        # Überprüfe, ob die Hauptvertices existieren
        additional_points = None
        if (x, y) in vert_dict and (nx, ny) in vert_dict:
            if (y == 0 or ny == 0):
                if (x < 0 or nx < 0):
                    additional_points = [(x-1, y), (x, y+1)]
                if (x > 0 or nx > 0):
                    additional_points = [(x+1, y), (x, y-1)]
            if (x == 0 or nx == 0):
                if (y < 0 or ny < 0):
                    additional_points = [(x+1, y), (x, y-1)]
                if (y > 0 or ny > 0):
                    additional_points = [(x-1, y), (x, y+1)]
            if (additional_points == None):
                additional_points = [(x-1, y), (x, y+1)]

            if additional_points is not None:
                for ap in additional_points:
                    if ap in vert_dict:
                        try:
                            # Erstelle das Dreieck
                            bm.faces.new([vert_dict[(x, y)], vert_dict[(nx, ny)], vert_dict[ap]])
                        except ValueError:
                            # Face existiert bereits
                            continue
    return bm


def add_diagonal_facesUL(bm, vert_dict, facesUL):
    for (x, y), (nx, ny) in facesUL:
        # Überprüfe, ob die Hauptvertices existieren
        additional_points = None
        if (x, y) in vert_dict and (nx, ny) in vert_dict:
            if (y == 0 or ny == 0):
                if (x < 0 or nx < 0):
                    # Bestimme die zusätzlichen Punkte für die Dreiecke
                    additional_points = [(x-1, y), (x, y-1)]
                if (x > 0 or nx > 0):
                    additional_points = [(x, y+1), (x+1, y)]
            if (y < 0 or ny < 0):
                if (x == 0 or nx == 0):
                    additional_points = [(x, y-1), (x-1, y)]
            if (y > 0 or ny > 0):
                if (x == 0 or nx == 0):
                    additional_points = [(x, y+1), (x+1, y)]
            if (additional_points == None):
                additional_points = [(x+1, y), (x, y+1)]
                
            if additional_points is not None:
                for ap in additional_points:
                    if ap in vert_dict:
                        try:
                            # Erstelle das Dreieck
                            bm.faces.new([vert_dict[(x, y)], vert_dict[(nx, ny)], vert_dict[ap]])
                        except ValueError:
                            # Face existiert bereits
                            continue
                        
    return bm


def copy_and_shift_mesh(bm, z_shift):
    new_verts = {}
    new_edges = {}

    original_verts = list(bm.verts)
    for vert in original_verts:
        new_vert = bm.verts.new((vert.co.x, vert.co.y, vert.co.z + z_shift))
        new_verts[vert] = new_vert

    print('Original vertices:', len(original_verts))
    print('New vertices:', len(new_verts))

    original_edges = list(bm.edges)
    for edge in original_edges:
        new_edge = bm.edges.new((new_verts[edge.verts[0]], new_verts[edge.verts[1]]))
        new_edges[edge] = new_edge

    print('Original edges:', len(original_edges))
    print('New edges:', len(new_edges))

    return bm

def find_boundary_vertices(bm):
    bm.verts.index_update()
    boundary_verts = set()
    for edge in bm.edges:
        if len(edge.link_faces) < 2:
            boundary_verts.add(edge.verts[0])
            boundary_verts.add(edge.verts[1])

    bm.verts.index_update()

    return list(boundary_verts)

def connect_meshes(bm, new_verts, new_edges):
    boundary_verts_original = find_boundary_vertices(bm)

    for v1 in boundary_verts_original:
        if v1 not in new_verts:
            print("Vertex nicht im new_verts-Dictionary gefunden:", v1)
            continue  # Überspringe diesen Vertex

        v2 = new_verts[v1]
        bm.edges.new([v1, v2])

    return bm


def setMaterial(obj):
    # Erstelle ein neues Material
    mat = bpy.data.materials.new(name="Benutzerdefiniertes_Material")

    # Aktiviere die Verwendung von Nodes (notwendig für die meisten Anpassungen)
    mat.use_nodes = True

    # Erhalte den Principled BSDF Shader Node
    principled_bsdf = mat.node_tree.nodes.get('Principled BSDF')
    if principled_bsdf is not None:
        # Setze die Basisfarbe (RGB)
        principled_bsdf.inputs['Base Color'].default_value = (1.0, 0.0, 0.0, 1)  # Rot

        # Setze die Rauheit
        principled_bsdf.inputs['Roughness'].default_value = 0.5  # Mittlere Rauheit

    # Weise das Material dem aktiven Objekt zu
    obj = bpy.context.object
    if obj and obj.type == 'MESH':
        if obj.data.materials:
            obj.data.materials[0] = mat
        else:
            obj.data.materials.append(mat)
            
    return obj

def animateMeshes(points_list, obj):
    start_frame = 1
    end_frame = len(points_list)
    bpy.context.scene.frame_start = start_frame
    bpy.context.scene.frame_end = end_frame

    # Stelle sicher, dass Shape Keys existieren
    if obj.data.shape_keys is None:
        obj.shape_key_add(name='Basis', from_mix=False)

    # Erstelle Shape Keys für alle Frames
    for frame_number, points in enumerate(points_list, start=start_frame):
        shape_key = obj.shape_key_add(name=f'Frame_{frame_number}', from_mix=False)
        for (vert, new_pos) in zip(shape_key.data, points):
            vert.co = new_pos

    # Animationsloop
    for frame_number in range(start_frame, end_frame + 1):
        bpy.context.scene.frame_set(frame_number)

        # Aktiviere den Shape Key für den aktuellen Frame und deaktiviere alle anderen
        for key_block in obj.data.shape_keys.key_blocks:
            if key_block.name == f'Frame_{frame_number}':
                key_block.value = 1.0
                key_block.keyframe_insert(data_path="value", frame=frame_number)
            else:
                key_block.value = 0.0
                key_block.keyframe_insert(data_path="value", frame=frame_number)


    # Optional: Zurück zum Start-Frame
    #bpy.context.scene.frame_set(start_frame)
    
def createFaces(bm):
    errorCounter = 0
    for face_verts in bm.verts:
        try:
            # Versuche, ein Face zu erstellen
            bm.faces.new(face_verts)
        except ValueError:
            errorCounter += 1
            continue

        print(errorCounter, ' Fehler')

    return bm

def createOuterFaces(bm):
    
    vert_dict = {v.index: v for v in bm.verts}

    # Finde offene Kanten, deren Länge größer als 1.0 ist
    offene_kanten = [kante for kante in bm.edges if len(kante.link_faces) < 2]
    missingEdges = [kante for kante in offene_kanten if kante.calc_length() > 1.0]
    
    for edge in missingEdges:
        p1, p2 = edge.verts

        # Finde Vertices, die mit p1 und p2 über eine andere Kante verbunden sind
        linkedVertsP1 = set(edge.other_vert(p1) for edge in p1.link_edges)
        linkedVertsP2 = set(edge.other_vert(p2) for edge in p2.link_edges)
        
        # Berechne die Schnittmenge der verbundenen Vertices
        gemeinsame_vertices = linkedVertsP1.intersection(linkedVertsP2)

        for v in gemeinsame_vertices:
            try:
                # Erstelle das Dreieck
                bm.faces.new([vert_dict[v.index], vert_dict[p1.index], vert_dict[p2.index]])
            except ValueError:
                # Face existiert bereits
                continue

def main():
    filename = "/home/alex/Dokumente/VsCode/Skripte/Python_Skript/Laser_Sound_Simulation/animation_data.json"
    points = import_data_from_json(filename)
    object_name = 'wavesurface'
    z_shift = 0.1

    # Erstelle ein neues Mesh und Objekt
    mesh = bpy.data.meshes.new(name="Mesh")
    obj = bpy.data.objects.new(object_name, mesh)
    bpy.context.collection.objects.link(obj)
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)
    bm = bmesh.new()

    bm = create_grid_mesh(points[0], bm, 20)
    print('create grid mesh finished')
    
    bm.to_mesh(mesh)
    mesh.update()
    
    bm.edges.ensure_lookup_table()
    bm.verts.ensure_lookup_table()
    
    createOuterFaces(bm)
    
    bm.to_mesh(mesh)
    mesh.update()
    
    bm.free()
    
    
    obj = setMaterial(obj)
    animateMeshes(points, obj)

if __name__ == "__main__":
    main()

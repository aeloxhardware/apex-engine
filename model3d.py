"""Modèle de données représentant un objet STL positionné sur le plateau d'impression."""
import os
import math
import struct
import uuid


class Model3D:
    """Représente un modèle STL individuel positionné sur le plateau."""
    def __init__(self, filepath, id_str=None):
        self.id = id_str if id_str else str(uuid.uuid4())
        self.filename = os.path.basename(filepath)
        self.filepath = filepath
        self.base_triangles = []
        
        # Vecteurs de transformation
        self.pos = [0.0, 0.0, 0.0]
        self.rot = [0.0, 0.0, 0.0]
        self.scale = [1.0, 1.0, 1.0]
        
        self.selected = False
        self.preview_triangles = None # Stocke l'aperçu dynamique de simplification
        
        # Caches pour l'accélération du rendu 3D
        self._intensity_cache_id = None
        self._triangle_intensities = []
        self._projected_polygons = [] # Requis pour la détection précise des clics (Hit Test)
        
        self.load_stl()

    def load_stl(self):
        """Lit et parse le fichier STL (ASCII ou Binaire)."""
        try:
            with open(self.filepath, 'rb') as f:
                header = f.read(80)
                if len(header) < 80: return
                is_ascii = False
                try:
                    if header.decode('utf-8', errors='ignore').strip().startswith('solid'):
                        is_ascii = True
                except Exception: pass
                
                if is_ascii:
                    f.seek(0)
                    content = f.read().decode('utf-8', errors='ignore')
                    current_tri = []
                    for line in content.splitlines():
                        line = line.strip()
                        if line.startswith('vertex'):
                            parts = line.split()
                            if len(parts) >= 4:
                                current_tri.append([float(parts[1]), float(parts[2]), float(parts[3])])
                        elif line.startswith('endloop'):
                            if len(current_tri) == 3:
                                self.base_triangles.append(current_tri)
                            current_tri = []
                else:
                    f.seek(80)
                    count_bytes = f.read(4)
                    if len(count_bytes) < 4: return
                    count = struct.unpack('<I', count_bytes)[0]
                    # Sous-échantillonnage dynamique si le modèle est trop lourd pour le moteur QPainter
                    step = int(count / 20000) + 1 if count > 20000 else 1
                    for i in range(count):
                        data = f.read(50)
                        if len(data) < 50: break
                        if i % step != 0: continue
                        t = struct.unpack('<fff fff fff fff H', data)
                        self.base_triangles.append([[t[3], t[4], t[5]], [t[6], t[7], t[8]], [t[9], t[10], t[11]]])
            
            self.center_origin()
        except Exception as e:
            print(f"Erreur chargement STL {self.filename}: {e}")

    def center_origin(self):
        """Centre le maillage sur son origine Z=0."""
        if not self.base_triangles: return
        xs = [v[0] for tri in self.base_triangles for v in tri]
        ys = [v[1] for tri in self.base_triangles for v in tri]
        zs = [v[2] for tri in self.base_triangles for v in tri]
        cx = (min(xs) + max(xs)) / 2.0
        cy = (min(ys) + max(ys)) / 2.0
        cz = min(zs)
        
        centered = []
        for tri in self.base_triangles:
            centered.append([[v[0]-cx, v[1]-cy, v[2]-cz] for v in tri])
        self.base_triangles = centered

    def get_world_triangles(self):
        """Applique l'échelle, la rotation et la translation aux triangles du modèle."""
        world_tris = []
        cx = math.cos(math.radians(self.rot[0]))
        sx = math.sin(math.radians(self.rot[0]))
        cy = math.cos(math.radians(self.rot[1]))
        sy = math.sin(math.radians(self.rot[1]))
        cz = math.cos(math.radians(self.rot[2]))
        sz = math.sin(math.radians(self.rot[2]))
        
        # Choisir le maillage d'aperçu de simplification ou le maillage d'origine
        source_triangles = self.preview_triangles if self.preview_triangles is not None else self.base_triangles
        
        for tri in source_triangles:
            new_tri = []
            for v in tri:
                # Échelle (Scale)
                x, y, z = v[0]*self.scale[0], v[1]*self.scale[1], v[2]*self.scale[2]
                # Rotation X
                y1 = y*cx - z*sx
                z1 = y*sx + z*cx
                # Rotation Y
                x2 = x*cy + z1*sy
                z2 = -x*sy + z1*cy
                # Rotation Z
                x3 = x2*cz - y1*sz
                y3 = x2*sz + y1*cz
                # Translation
                new_tri.append([x3+self.pos[0], y3+self.pos[1], z2+self.pos[2]])
            world_tris.append(new_tri)
        return world_tris

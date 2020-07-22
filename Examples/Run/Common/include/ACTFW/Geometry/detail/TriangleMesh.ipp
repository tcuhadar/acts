// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename T>
void TriangleMesh<T>::surface(const Surface& surface, const GeometryContext& gctx) {
  const Vector3D position(0.0, 0.0, 0.0);
  const Vector3D direction(0.3, 0.4, 0.0);
  const BoundaryCheck& bcheck(true);

  std::cout << "TriangleMesh::surface " << surface.type() << std::endl;
  m_surface_intersection = surface.intersect(gctx, position, direction, bcheck);
  std::cout << "  " << (bool)m_surface_intersection << std::endl;
  m_surface_count++;
  m_surface_type = surface.type();
  if (m_surface_intersection) {
    m_surface_intersection_count++;
  }
}

template <typename T>
void TriangleMesh<T>::vertex(const Vector3D& /* vtx */, ColorRGB /* color */) {
  std::cout << "TriangleMesh::vertex() NOT IMPLEMENTED" << std::endl;
  // m_vertexColors[m_vertices.size()] = color;
  // m_vertices.push_back(vtx.template cast<ValueType>());
}

template <typename T>
void TriangleMesh<T>::line(const Vector3D& /* a */, const Vector3D& /* b */,
                               ColorRGB /* color */) {
  std::cout << "TriangleMesh::line() NOT IMPLEMENTED" << std::endl;
  // if (color != ColorRGB{0, 0, 0}) {
  //   m_lineColors[m_lines.size()] = color;
  // }
  // // not implemented
  // vertex(a, color);
  // vertex(b, color);
  // m_lines.push_back({m_vertices.size() - 2, m_vertices.size() - 1});
}

template <typename T>
void TriangleMesh<T>::face(const std::vector<Vector3D>& /* vtxs */,
                               ColorRGB /* color */) {
  std::cout << "TriangleMesh::face() NOT IMPLEMENTED" << std::endl;
  // if (color != ColorRGB{0, 0, 0}) {
  //   m_faceColors[m_faces.size()] = color;
  // }
  // FaceType idxs;
  // idxs.reserve(vtxs.size());
  // for (const auto& vtx : vtxs) {
  //   vertex(vtx, color);
  //   idxs.push_back(m_vertices.size() - 1);
  // }
  // m_faces.push_back(std::move(idxs));
}

template <typename T>
void TriangleMesh<T>::faces(const std::vector<Vector3D>& vtxs,
                                const std::vector<FaceType>& faces,
                                ColorRGB color) {
  // std::cout << "TriangleMesh::faces()" << std::endl;
  // No faces given - call the face() method
  if (faces.empty()) {
    face(vtxs, color);
  } else {
    // if (color != ColorRGB{0, 0, 0}) {
    //   m_faceColors[m_faces.size()] = color;
    // }
    auto vtxoffs = m_vertices.size();
    // if (color != ColorRGB{0, 0, 0}) {
    //   m_vertexColors[m_vertices.size()] = color;
    // }
    m_vertices.insert(m_vertices.end(), vtxs.begin(), vtxs.end());
    for (const auto& face : faces) {
      if (face.size() == 2) {
        std::cout << "TriangleMesh::faces() - line NOT IMPLEMENTED" << std::endl;
        // m_lines.push_back({face[0] + vtxoffs, face[2] + vtxoffs});
      } else if (face.size() > 3) {
        std::cout << "TriangleMesh::faces() - not a triangle " << face.size() << std::endl;
      } else {
        FaceType rawFace = face;
        std::transform(rawFace.begin(), rawFace.end(), rawFace.begin(),
                       [&](size_t& iv) { return (iv + vtxoffs); });
        m_faces.push_back(rawFace);
        m_surface_intersect[rawFace] = m_surface_intersection;
        m_surface_intersect_count[rawFace] = m_surface_count;
        m_surface_intersect_type[rawFace] = m_surface_type;

        const Vector3D position(0.0, 0.0, 0.0);
        const Vector3D direction(0.3, 0.4, 0.0);
        m_face_intersect[rawFace] = intersect(position, direction, vtxs[face[0]], vtxs[face[1]], vtxs[face[2]]);
        if (m_face_intersect[rawFace]) {
          m_face_intersection_count++;
        }
      }
    }
  }
}

template <typename T>
void TriangleMesh<T>::write(const std::string& /* path */) const {
  // std::cout << "TriangleMesh::write(" << path << ") NOT IMPLEMENTED" << std::endl;
  // std::ofstream os;
  // std::string objectpath = path;
  // if (not IVisualization::hasExtension(objectpath)) {
  //   objectpath += std::string(".obj");
  // }
  // os.open(objectpath);
  // std::string mtlpath = objectpath;
  // IVisualization::replaceExtension(mtlpath, ".mtl");
  // os << "mtllib " << mtlpath << "\n";
  // std::ofstream mtlos;
  // mtlos.open(mtlpath);
  // write(os, mtlos);
  // os.close();
  // mtlos.close();
}

template <typename T>
void TriangleMesh<T>::write(std::ostream& /* os */) const {
  std::cout << "TriangleMesh " << m_vertices.size() << " vertices, "
                               << m_faces.size() << " faces, "
                               << m_surface_count << " surfaces, "
                               << m_surface_intersection_count << " surface intersects, "
                               << m_face_intersection_count << " face intersects"
                               << std::endl;

  std::ofstream tos;
  tos.open("TriangleMesh.obj");

  for (const VertexType& vtx : m_vertices) {
    tos << "v " << std::setprecision(m_outputPrecision)
        << m_outputScalor * vtx.x() << " " << m_outputScalor * vtx.y() << " "
        << m_outputScalor * vtx.z() << "\n";
  }

  // Surface
  // for (const FaceType fc : m_faces) {
  //   SurfaceIntersection si = m_surface_intersect.at(fc);
  //   if (si) {
  //     tos << "f";
  //     for (size_t i = 0; i < fc.size(); i++) {
  //       tos << " " << fc[i] + 1;
  //     }
  //     tos << "  # " << m_surface_intersect_count.at(fc);
  //     tos << "\n";
  //   }
  // }

  // Faces
  for (const FaceType fc : m_faces) {
    bool fi = (bool)m_face_intersect.at(fc);
    unsigned int type = m_surface_intersect_type.at(fc);
    if (fi) {
      tos << "# " << m_surface_intersect_count.at(fc) << " "
          << type << " "
          << (bool)m_surface_intersect.at(fc)
          << " " << (bool)m_face_intersect.at(fc)
          << "\n";
      tos << "f";
      for (size_t i = 0; i < fc.size(); i++) {
        tos << " " << fc[i] + 1;
      }
      tos << "\n";
    }
  }


  tos.close();
  // std::stringstream sterile;
  // write(os, sterile);
}

template <typename T>
void TriangleMesh<T>::write(std::ostream& /* os */, std::ostream& /* mos */) const {
  // std::cout << "TriangleMesh::write() NOT IMPLEMENTED" << std::endl;
  // std::map<std::string, bool> materials;

  // auto mixColor = [&](const ColorRGB& color) -> std::string {
  //   std::string materialName;
  //   materialName = "material_";
  //   materialName += std::to_string(color[0]) + std::string("_");
  //   materialName += std::to_string(color[1]) + std::string("_");
  //   materialName += std::to_string(color[2]);

  //   if (materials.find(materialName) == materials.end()) {
  //     mos << "newmtl " << materialName << "\n";
  //     std::vector<std::string> shadings = {"Ka", "Kd", "Ks"};
  //     for (const auto& shd : shadings) {
  //       mos << shd << " " << std::to_string(color[0] / 256.) << " ";
  //       mos << std::to_string(color[1] / 256.) << " ";
  //       mos << std::to_string(color[2] / 256.) << " "
  //           << "\n";
  //     }
  //     mos << "\n";
  //   }
  //   return std::string("usemtl ") + materialName;
  // };

  // size_t iv = 0;
  // ColorRGB lastVertexColor = {0, 0, 0};
  // for (const VertexType& vtx : m_vertices) {
  //   if (m_vertexColors.find(iv) != m_vertexColors.end()) {
  //     auto color = m_vertexColors.find(iv)->second;
  //     if (color != lastVertexColor) {
  //       os << mixColor(color) << "\n";
  //       lastVertexColor = color;
  //     }
  //   }

  //   os << "v " << std::setprecision(m_outputPrecision)
  //      << m_outputScalor * vtx.x() << " " << m_outputScalor * vtx.y() << " "
  //      << m_outputScalor * vtx.z() << "\n";
  //   ++iv;
  // }
  // size_t il = 0;
  // ColorRGB lastLineColor = {0, 0, 0};
  // for (const LineType& ln : m_lines) {
  //   if (m_lineColors.find(il) != m_lineColors.end()) {
  //     auto color = m_lineColors.find(il)->second;
  //     if (color != lastLineColor) {
  //       os << mixColor(color) << "\n";
  //       lastLineColor = color;
  //     }
  //   }
  //   os << "l " << ln.first + 1 << " " << ln.second + 1 << "\n";
  //   ++il;
  // }
  // size_t is = 0;
  // ColorRGB lastFaceColor = {0, 0, 0};
  // for (const FaceType& fc : m_faces) {
  //   if (m_faceColors.find(is) != m_faceColors.end()) {
  //     auto color = m_faceColors.find(is)->second;
  //     if (color != lastFaceColor) {
  //       os << mixColor(color) << "\n";
  //       lastFaceColor = color;
  //     }
  //   }
  //   os << "f";
  //   for (size_t i = 0; i < fc.size(); i++) {
  //     os << " " << fc[i] + 1;
  //   }
  //   os << "\n";
  //   ++is;
  // }
}

template <typename T>
void TriangleMesh<T>::clear() {
  // std::cout << "TriangleMesh::clear() NOT IMPLEMENTED" << std::endl;
  // m_vertices.clear();
  // m_faces.clear();
  // m_lines.clear();
  // m_lineColors.clear();
  // m_vertexColors.clear();
  // m_faceColors.clear();
}

// from: https://stackoverflow.com/questions/42740765/intersection-between-line-and-triangle-in-3d
template <typename T>
bool TriangleMesh<T>::intersect(const Vector3D& orig,
                                const Vector3D& direction,
                                const Vector3D& v0,
                                const Vector3D& v1,
                                const Vector3D& v2) const {
    Vector3D v0v1 = v1 - v0; 
    Vector3D v0v2 = v2 - v0; 
    Vector3D dir = direction.normalized();
    Vector3D pvec = dir.cross(v0v2); 
    float det = v0v1.dot(pvec); 
    // ray and triangle are parallel if det is close to 0
    if (fabs(det) < 1e-6) return false; 
    double invDet = 1 / det; 
 
    Vector3D tvec = orig - v0; 
    double u = tvec.dot(pvec) * invDet; 
    if (u < 0 || u > 1) return false; 
 
    Vector3D qvec = tvec.cross(v0v1); 
    double v = dir.dot(qvec) * invDet; 
    if (v < 0 || u + v > 1) return false; 
 
    double t = v0v2.dot(qvec) * invDet; 
    if (t < 0) return false;
    std::cout << "det = "<< det << " t =  " << t << " u = " << u << " v = " << v << " u+v = " << (u+v) << std::endl;
 
    return true; 
}

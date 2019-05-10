// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename entity_t, typename value_t, size_t DIM>
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::AxisAlignedBoundingBox(
    const entity_t* entity, const vertex_type& vmin, const vertex_type& vmax)
    : m_entity(entity),
      m_vmin(vmin),
      m_vmax(vmax),
      m_center((vmin + vmax) / 2.),
      m_width(vmax - vmin),
      m_iwidth(1 / m_width) {}

template <typename entity_t, typename value_t, size_t DIM>
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::AxisAlignedBoundingBox(
    const entity_t* entity, const vertex_type& center, const Size& size)
    : m_entity(entity),
      m_vmin(center - size.get() * 0.5),
      m_vmax(center + size.get() * 0.5),
      m_center(center),
      m_width(size.get()),
      m_iwidth(1 / m_width) {}

template <typename entity_t, typename value_t, size_t DIM>
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::AxisAlignedBoundingBox(
    const std::vector<self_t*>& boxes, vertex_array_type envelope)
    : m_entity(nullptr) {
  assert(boxes.size() > 1);

  for (size_t i = 0; i < boxes.size(); i++) {
    if (i < boxes.size() - 1) {
      // set next on i to i+1
      boxes[i]->setSkip(boxes[i + 1]);
    } else {
      // make sure last is set to nullptr, this marks end
      // boxes[i]->m_next = nullptr;
      boxes[i]->setSkip(nullptr);
    }
  }

  m_left_child = boxes.front();
  m_right_child = boxes.back();
  m_skip = nullptr;

  std::tie(m_vmin, m_vmax) = wrap(boxes, envelope);

  m_center = (m_vmin + m_vmax) / 2.;
  m_width = m_vmax - m_vmin;
  m_iwidth = 1 / m_width;
}

template <typename entity_t, typename value_t, size_t DIM>
std::pair<
    typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::vertex_type,
    typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::vertex_type>
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::wrap(
    const std::vector<const self_t*>& boxes, vertex_array_type envelope) {
  assert(boxes.size() > 1);
  // figure out extent of boxes
  // use array for Eigen coefficient wise min/max
  vertex_array_type vmax(
      vertex_array_type::Constant(std::numeric_limits<value_type>::lowest()));
  vertex_array_type vmin(
      vertex_array_type::Constant(std::numeric_limits<value_type>::max()));

  for (size_t i = 0; i < boxes.size(); i++) {
    vmin = vmin.min(boxes[i]->min().array());
    vmax = vmax.max(boxes[i]->max().array());
  }

  vmax += envelope;
  vmin -= envelope;

  return {vmin, vmax};
}

template <typename entity_t, typename value_t, size_t DIM>
std::pair<
    typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::vertex_type,
    typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::vertex_type>
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::wrap(
    const std::vector<self_t*>& boxes, vertex_array_type envelope) {
  assert(boxes.size() > 1);
  std::vector<const self_t*> box_ptrs;
  box_ptrs.reserve(boxes.size());
  std::transform(boxes.begin(), boxes.end(), std::back_inserter(box_ptrs),
                 [](const auto* box) { return box; });
  return wrap(box_ptrs, envelope);
}

template <typename entity_t, typename value_t, size_t DIM>
std::pair<
    typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::vertex_type,
    typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::vertex_type>
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::wrap(
    const std::vector<self_t>& boxes, vertex_array_type envelope) {
  assert(boxes.size() > 1);
  std::vector<const self_t*> box_ptrs;
  box_ptrs.reserve(boxes.size());
  std::transform(boxes.begin(), boxes.end(), std::back_inserter(box_ptrs),
                 [](auto& box) { return &box; });
  return wrap(box_ptrs, envelope);
}

template <typename entity_t, typename value_t, size_t DIM>
bool Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::intersect(
    const vertex_type& point) const {
  vertex_array_type t = (point - m_vmin).array() * m_iwidth;
  return t.minCoeff() >= 0 && t.maxCoeff() < 1;
}

template <typename entity_t, typename value_t, size_t DIM>
void Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::setSkip(
    self_t* skip) {
  // set next on this
  m_skip = skip;
  // find last child and set its skip
  if (m_right_child != nullptr) {
    m_right_child->setSkip(skip);
  }
}

template <typename entity_t, typename value_t, size_t DIM>
const Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>*
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::getLeftChild() const {
  return m_left_child;
}

template <typename entity_t, typename value_t, size_t DIM>
const Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>*
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::getSkip() const {
  return m_skip;
}

template <typename entity_t, typename value_t, size_t DIM>
bool Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::hasEntity() const {
  return m_entity != nullptr;
}

template <typename entity_t, typename value_t, size_t DIM>
const entity_t* Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::entity()
    const {
  return m_entity;
}

template <typename entity_t, typename value_t, size_t DIM>
void Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::setEntity(
    const entity_t* entity) {
  m_entity = entity;
}

template <typename entity_t, typename value_t, size_t DIM>
const typename Acts::AxisAlignedBoundingBox<entity_t, value_t,
                                            DIM>::vertex_type&
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::center() const {
  return m_center;
}

template <typename entity_t, typename value_t, size_t DIM>
const typename Acts::AxisAlignedBoundingBox<entity_t, value_t,
                                            DIM>::vertex_type&
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::min() const {
  return m_vmin;
}

template <typename entity_t, typename value_t, size_t DIM>
const typename Acts::AxisAlignedBoundingBox<entity_t, value_t,
                                            DIM>::vertex_type&
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::max() const {
  return m_vmax;
}

template <typename entity_t, typename value_t, size_t DIM>
std::ostream& Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::toStream(
    std::ostream& os) const {
  os << "AABB(ctr=(";

  for (size_t i = 0; i < DIM; i++) {
    if (i > 0) {
      os << ", ";
    }
    os << m_center[i];
  }

  os << ") vmin=(";
  for (size_t i = 0; i < DIM; i++) {
    if (i > 0) {
      os << ", ";
    }
    os << m_vmin[i];
  }

  os << ") vmax=(";

  for (size_t i = 0; i < DIM; i++) {
    if (i > 0) {
      os << ", ";
    }
    os << m_vmax[i];
  }

  os << "))";

  return os;
}

template <typename entity_t, typename value_t, size_t DIM>
template <size_t D, std::enable_if_t<D == 3, int>>
std::pair<
    typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::vertex_type,
    typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::vertex_type>
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::transformVertices(
    const transform_type& trf) const {
  // we need to enumerate all the vertices, transform,
  // and then recalculate min and max

  std::array<vertex_type, 8> vertices({{
      {m_vmin.x(), m_vmin.y(), m_vmin.z()},
      {m_vmin.x(), m_vmax.y(), m_vmin.z()},
      {m_vmax.x(), m_vmax.y(), m_vmin.z()},
      {m_vmax.x(), m_vmin.y(), m_vmin.z()},
      {m_vmin.x(), m_vmin.y(), m_vmax.z()},
      {m_vmin.x(), m_vmax.y(), m_vmax.z()},
      {m_vmax.x(), m_vmax.y(), m_vmax.z()},
      {m_vmax.x(), m_vmin.y(), m_vmax.z()},
  }});

  vertex_type vmin = trf * vertices[0];
  vertex_type vmax = trf * vertices[0];

  for (size_t i = 1; i < 8; i++) {
    const vertex_type vtx = trf * vertices[i];
    vmin = vmin.cwiseMin(vtx);
    vmax = vmax.cwiseMax(vtx);
  }

  return {vmin, vmax};
}

template <typename entity_t, typename value_t, size_t DIM>
template <size_t D, std::enable_if_t<D == 2, int>>
std::pair<
    typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::vertex_type,
    typename Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::vertex_type>
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::transformVertices(
    const transform_type& trf) const {
  // we need to enumerate all the vertices, transform,
  // and then recalculate min and max

  std::array<vertex_type, 4> vertices({{{m_vmin.x(), m_vmin.y()},
                                        {m_vmin.x(), m_vmax.y()},
                                        {m_vmax.x(), m_vmax.y()},
                                        {m_vmax.x(), m_vmin.y()}}});

  vertex_type vmin = trf * vertices[0];
  vertex_type vmax = trf * vertices[0];

  for (size_t i = 1; i < 4; i++) {
    const vertex_type vtx = trf * vertices[i];
    vmin = vmin.cwiseMin(vtx);
    vmax = vmax.cwiseMax(vtx);
  }

  return {vmin, vmax};
}

template <typename entity_t, typename value_t, size_t DIM>
void Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::transform(
    const transform_type& trf) {
  std::tie(m_vmin, m_vmax) = transformVertices(trf);
}

template <typename entity_t, typename value_t, size_t DIM>
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>
Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::transformed(
    const transform_type& trf) const {
  vertex_type vmin, vmax;
  std::tie(vmin, vmax) = transformVertices(trf);
  return self_t(m_entity, vmin, vmax);
}

template <typename entity_t, typename value_t, size_t DIM>
template <size_t D, std::enable_if_t<D == 3, int>>
void Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::draw(
    IVisualization& helper, std::array<int, 3> color,
    const transform_type& trf) const {
  static_assert(DIM == 3, "PLY output only supported in 3D");

  const vertex_type& vmin = m_vmin;
  const vertex_type& vmax = m_vmax;

  auto write = [&](const vertex_type& a, const vertex_type& b,
                   const vertex_type& c, const vertex_type& d) {
    helper.face(std::vector<vertex_type>({trf * a, trf * b, trf * c, trf * d}),
                color);
  };

  write({vmin.x(), vmin.y(), vmin.z()}, {vmin.x(), vmax.y(), vmin.z()},
        {vmin.x(), vmax.y(), vmax.z()}, {vmin.x(), vmin.y(), vmax.z()});

  write({vmax.x(), vmin.y(), vmin.z()}, {vmax.x(), vmax.y(), vmin.z()},
        {vmax.x(), vmax.y(), vmax.z()}, {vmax.x(), vmin.y(), vmax.z()});

  write({vmin.x(), vmin.y(), vmin.z()}, {vmax.x(), vmin.y(), vmin.z()},
        {vmax.x(), vmin.y(), vmax.z()}, {vmin.x(), vmin.y(), vmax.z()});

  write({vmin.x(), vmax.y(), vmin.z()}, {vmax.x(), vmax.y(), vmin.z()},
        {vmax.x(), vmax.y(), vmax.z()}, {vmin.x(), vmax.y(), vmax.z()});

  write({vmin.x(), vmin.y(), vmin.z()}, {vmax.x(), vmin.y(), vmin.z()},
        {vmax.x(), vmax.y(), vmin.z()}, {vmin.x(), vmax.y(), vmin.z()});

  write({vmin.x(), vmin.y(), vmax.z()}, {vmax.x(), vmin.y(), vmax.z()},
        {vmax.x(), vmax.y(), vmax.z()}, {vmin.x(), vmax.y(), vmax.z()});
}

template <typename entity_t, typename value_t, size_t DIM>
template <size_t D, std::enable_if_t<D == 2, int>>
std::ostream& Acts::AxisAlignedBoundingBox<entity_t, value_t, DIM>::svg(
    std::ostream& os, value_type w, value_type h, value_type unit,
    std::string label, std::string fillcolor) const {
  static_assert(DIM == 2, "SVG is only supported in 2D");

  vertex_type mid(w / 2., h / 2.);

  using transform_t = Eigen::Transform<value_t, DIM, Eigen::Affine>;

  transform_t trf = transform_t::Identity();
  trf.translate(mid);
  trf = trf * Eigen::Scaling(vertex_type(1, -1));
  trf.scale(unit);

  auto draw_point = [&](const vertex_type& p_, std::string color, size_t r) {
    vertex_type p = trf * p_;
    os << "<circle ";
    os << "cx=\"" << p.x() << "\" cy=\"" << p.y() << "\" r=\"" << r << "\"";
    os << " fill=\"" << color << "\"";
    os << "/>\n";
  };

  auto draw_rect = [&](const vertex_type& center_, const vertex_type& size_,
                       std::string color) {
    vertex_type size = size_ * unit;
    vertex_type center = trf * center_ - size * 0.5;

    os << "<rect ";
    os << "x=\"" << center.x() << "\" y=\"" << center.y() << "\" ";
    os << "width=\"" << size.x() << "\" height=\"" << size.y() << "\"";
    os << " fill=\"" << color << "\"";
    os << "/>\n";
  };

  auto draw_text = [&](const vertex_type& center_, std::string text,
                       std::string color, size_t size) {
    vertex_type center = trf * center_;
    os << "<text dominant-baseline=\"middle\" text-anchor=\"middle\" ";
    os << "fill=\"" << color << "\" font-size=\"" << size << "\" ";
    os << "x=\"" << center.x() << "\" y=\"" << center.y() << "\">";
    os << text << "</text>\n";
  };

  draw_rect(m_center, m_width, fillcolor);
  draw_point(m_vmin, "black", 2);
  draw_point(m_vmax, "black", 2);
  draw_text(m_center, label, "white", 10);

  return os;
}

template <typename box_t>
box_t* octree_inner(std::vector<std::unique_ptr<box_t>>& store,
                    size_t max_depth,
                    typename box_t::vertex_array_type envelope,
                    const std::vector<box_t*>& lprims, size_t depth) {
  using vertex_type = typename box_t::vertex_type;

  assert(lprims.size() > 0);
  if (lprims.size() == 1) {
    // just return
    return lprims.front();
  }

  if (depth >= max_depth) {
    // just wrap them all up
    auto bb = std::make_unique<box_t>(lprims, envelope);
    store.push_back(std::move(bb));
    return store.back().get();
  }

  std::array<std::vector<box_t*>, 8> octants;
  // calc center of boxes
  vertex_type vmin, vmax;
  std::tie(vmin, vmax) = box_t::wrap(lprims);
  vertex_type glob_ctr = (vmin + vmax) / 2.;

  for (auto* box : lprims) {
    vertex_type ctr = box->center() - glob_ctr;
    if (ctr.x() < 0 && ctr.y() < 0 && ctr.z() < 0) {
      octants[0].push_back(box);
      continue;
    }
    if (ctr.x() > 0 && ctr.y() < 0 && ctr.z() < 0) {
      octants[1].push_back(box);
      continue;
    }
    if (ctr.x() < 0 && ctr.y() > 0 && ctr.z() < 0) {
      octants[2].push_back(box);
      continue;
    }
    if (ctr.x() > 0 && ctr.y() > 0 && ctr.z() < 0) {
      octants[3].push_back(box);
      continue;
    }

    if (ctr.x() < 0 && ctr.y() < 0 && ctr.z() > 0) {
      octants[4].push_back(box);
      continue;
    }
    if (ctr.x() > 0 && ctr.y() < 0 && ctr.z() > 0) {
      octants[5].push_back(box);
      continue;
    }
    if (ctr.x() < 0 && ctr.y() > 0 && ctr.z() > 0) {
      octants[6].push_back(box);
      continue;
    }
    if (ctr.x() > 0 && ctr.y() > 0 && ctr.z() > 0) {
      octants[7].push_back(box);
      continue;
    }

    // not in any quadrant (numerics probably)
    octants[0].push_back(box);
  }

  std::vector<box_t*> sub_octs;
  for (const auto& sub_prims : octants) {
    if (sub_prims.size() <= 8) {
      if (sub_prims.size() < 1) {
        // done
      } else if (sub_prims.size() == 1) {
        sub_octs.push_back(sub_prims.front());
      } else {
        store.push_back(std::make_unique<box_t>(sub_prims, envelope));
        sub_octs.push_back(store.back().get());
      }
    } else {
      // recurse
      sub_octs.push_back(
          octree_inner(store, max_depth, envelope, sub_prims, depth + 1));
    }
  }

  if (sub_octs.size() == 1) {
    return sub_octs.front();
  }

  auto bb = std::make_unique<box_t>(sub_octs, envelope);
  store.push_back(std::move(bb));
  return store.back().get();
}

template <typename box_t>
box_t* Acts::make_octree(std::vector<std::unique_ptr<box_t>>& store,
                         const std::vector<box_t*>& prims, size_t max_depth,
                         typename box_t::value_type envelope1) {
  static_assert(box_t::dim == 3, "Octree can only be created in 3D");

  using vertex_array_type = typename box_t::vertex_array_type;

  vertex_array_type envelope(vertex_array_type::Constant(envelope1));

  box_t* top = octree_inner(store, max_depth, envelope, prims, 0);
  return top;
}

template <typename T, typename U, size_t V>
std::ostream& operator<<(std::ostream& os,
                         const Acts::AxisAlignedBoundingBox<T, U, V>& box) {
  box.dump(os);
  return os;
}
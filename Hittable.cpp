
#include "Hittable.h"

// Sphere
bool Sphere::Hit(const Ray& ray, HitRecord *hit_record) const {
    // TODO: Add your code here.

    float A ;
    float B ;
    float C ;
    float dis ;
    float t = 0;;
    float t_1 ;
    float t_2 ;
    float smaller_t;
    glm::vec3 normal;
    glm::vec3 position;
    glm::vec3 unit_norm;
    glm::vec3 unit_ray_d;
    
   
    A =  glm::dot(ray.d,ray.d);
    B =  2 * glm::dot(ray.o-o_,ray.d);
    C =  glm::dot(ray.o-o_,ray.o-o_) - pow(r_,2);
    dis = sqrt(pow(B,2)-4*C);
    

    
    if (dis >= 0){
        t_1 = (-B - dis)/2;
        t_2 = (-B + dis)/2; 


        if (t_1 > 0 && t_2 > 0){
            t = std::min(t_1,t_2);
        }
        else{
            if (t_1 > 0){
                t = t_1;
            }
            else if (t_2 > 0){
                t = t_2;
            }
            
        }

        if (t == 0){
            return false;
        }
        
        unit_ray_d = glm::normalize(ray.d);
        position = ray.At(t);
        normal =  position - o_;
        unit_norm = glm::normalize(normal);

        hit_record->material = material_;
        hit_record->in_direction = unit_ray_d;
        hit_record->normal = unit_norm;
        hit_record->position = position;
        hit_record->distance = glm::length( ray.d*t);
        
        hit_record->reflection = glm::normalize(unit_ray_d -2 * glm::dot(unit_norm,unit_ray_d) * unit_norm );
       
        return true;
    }
    return false;
    
    
} 


// Quadric
bool Quadric::Hit(const Ray& ray, HitRecord *hit_record) const {
    // TODO: Add your code here.
 
    glm::vec4 O =  glm::vec4(ray.o,  1);
    glm::vec4 D =  glm::vec4(ray.d,  0);
    glm::vec3 norm_normal;
    glm::vec3 norm_ray_d;
    float t_1 ;
    float t_2 ;
    //std::cout << "quadric" << std::endl;
    float a = glm::dot(D, A_ * D);
    float b = 2 * glm::dot(O, A_ * D);
    float c = glm::dot(O, A_ * O);
    float t = 0;
    
    glm::vec3 position;
    glm::vec3 normal;
    float discrim = pow(b,2) - 4 * a * c;

    
   
    if (discrim > 0){
        t_1 = (-b + sqrt(discrim))/(2*a);
        t_2 = (-b - sqrt(discrim))/(2*a);
        if (t_1 > 0 && t_2 >=0){
            t = std::min(t_1,t_2);
        }
        else{
            if (t_1 > 0){
                t = t_1;
            }
            else if (t_2 > 0){
                t = t_2;
            }
        }
    }

    if (t == 0){
        return false;
    }

    else if (discrim == 0){

        t = -b/(2*a);   
    }

    else {
        //std::cout<<"on9" <<std::endl;
 
        return false;
    }

    //std::cout<<"work" <<std::endl;
    position = ray.At(t);
    normal = (glm::transpose(A_) + A_) * (O+D *t);
    norm_normal = glm::normalize(normal);
    norm_ray_d = glm::normalize(ray.d);
    hit_record->material = material_;
    hit_record->in_direction = norm_ray_d;
    hit_record->position = position;
    hit_record->distance = glm::length(t * ray.d);
    hit_record->normal = norm_normal;
    hit_record->reflection = glm::normalize(norm_ray_d -2 * glm::dot(norm_normal,norm_ray_d) * norm_normal);

    return true;
   
}

// Triangle
bool Triangle::Hit(const Ray& ray, HitRecord *hit_record) const {
    // TODO: Add your code here.
    bool ret = false;
    //std::cout << "Triangle" << std::endl;
    glm::vec3 centroid;
    glm::vec3 normal;
    glm::vec3 position;

    glm::vec3 cross_1;
    glm::vec3 cross_2;
    glm::vec3 cross_3;

    glm::vec3 norm_normal;
    glm::vec3 norm_ray_d;

    float t;
    float alpha_1;
    float alpha_2;
    float alpha_3;
    float denominator;
    float triangle_area;
    float left_triangle_area;
    float right_triangle_area;
    float bottom_triangle_area;
    float D;

    glm::vec3 v_a_c_a = glm::cross((b_-a_),(c_-a_));

    centroid = (a_+b_+c_)/3.0f;
    
    D = glm::dot(v_a_c_a , a_);
    denominator = glm::dot(v_a_c_a,ray.d);
    if (denominator == 0){
        return false;
    }
    
    t = - (glm::dot(v_a_c_a,ray.o)-D)/denominator;
    if (t <= 0) {
        return false;
    }
    position =ray.At(t);

    
  
    cross_1 = (glm::cross((b_-a_),(position-a_)));
    cross_2 = (glm::cross((position-a_),(c_-a_)));
    cross_3 = (glm::cross((position-c_),(b_-c_)));
    
    triangle_area = glm::length(glm::cross((b_-a_),(c_-a_)))/2;
    left_triangle_area = glm::length(cross_1)/2;
    bottom_triangle_area = glm::length(cross_2)/2; 
    right_triangle_area = glm::length(cross_3)/2; 
    alpha_1 = left_triangle_area/triangle_area;
    alpha_2 = bottom_triangle_area/triangle_area;
    alpha_3 = right_triangle_area/triangle_area;

    if (phong_interpolation_){

        normal = alpha_1 * n_c_+ alpha_2 * n_b_+ alpha_3 * n_a_;
    }
    else{
        normal = v_a_c_a;
    }
    
    

    
      // check if inside 

    if (! (glm::dot(cross_1,normal) >= 0 && glm::dot(cross_2,normal) >= 0 && glm::dot(cross_3,normal) >= 0) ){
        return false;
    }
 

        
    norm_normal = glm::normalize(normal);
    norm_ray_d = glm::normalize(ray.d);
    hit_record->normal =norm_normal;
    hit_record->in_direction = norm_ray_d;
    hit_record->position = position;
    hit_record->distance = glm::length(t * ray.d);
    hit_record->reflection = glm::normalize(ray.d -2 * glm::dot(norm_normal,norm_ray_d) * norm_normal);
 

    
    return true;
}

// ---------------------------------------------------------------------------------------------
// ------------------------------ no need to change --------------------------------------------
// ---------------------------------------------------------------------------------------------

// CompleteTriangle
bool CompleteTriangle::Hit(const Ray& ray, HitRecord *hit_record) const {
    bool ret = triangle_.Hit(ray, hit_record);
    if (ret) {
        hit_record->material = material_;
    }
    return ret;
}


// Mesh
Mesh::Mesh(const std::string& file_path,
           const Material& material,
           bool phong_interpolation):
           ply_data_(file_path), material_(material), phong_interpolation_(phong_interpolation) {
    std::vector<std::array<double, 3>> v_pos = ply_data_.getVertexPositions();
    vertices_.resize(v_pos.size());

    for (int i = 0; i < vertices_.size(); i++) {
        vertices_[i] = Point(v_pos[i][0], v_pos[i][1], v_pos[i][2]);
    }

    f_ind_ = ply_data_.getFaceIndices();

    // Calc face normals
    for (const auto& face : f_ind_) {
        Vec normal = glm::normalize(glm::cross(vertices_[face[1]] - vertices_[face[0]], vertices_[face[2]] - vertices_[face[0]]));
        face_normals_.emplace_back(normal);
    }

    // Calc vertex normals
    vertex_normals_.resize(vertices_.size(), Vec(0.f, 0.f, 0.f));
    for (int i = 0; i < f_ind_.size(); i++) {
        for (int j = 0; j < 3; j++) {
            vertex_normals_[f_ind_[i][j]] += face_normals_[i];
        }
    }
    for (auto& vertex_normal : vertex_normals_) {
        vertex_normal = glm::normalize(vertex_normal);
    }

    // Construct hittable triangles
    for (const auto& face : f_ind_) {
        triangles_.emplace_back(vertices_[face[0]], vertices_[face[1]], vertices_[face[2]],
                                vertex_normals_[face[0]], vertex_normals_[face[1]], vertex_normals_[face[2]],
                                phong_interpolation_);
    }

    // Calc bounding box
    Point bbox_min( 1e5f,  1e5f,  1e5f);
    Point bbox_max(-1e5f, -1e5f, -1e5f);
    for (const auto& vertex : vertices_) {
        bbox_min = glm::min(bbox_min, vertex - 1e-3f);
        bbox_max = glm::max(bbox_max, vertex + 1e-3f);
    }

    // Build Octree
    tree_nodes_.emplace_back(new OctreeNode());
    tree_nodes_.front()->bbox_min = bbox_min;
    tree_nodes_.front()->bbox_max = bbox_max;

    root_ = tree_nodes_.front().get();
    for (int i = 0; i < f_ind_.size(); i++) {
        InsertFace(root_, i);
    }
}

bool Mesh::Hit(const Ray& ray, HitRecord *hit_record) const {
    const bool brute_force = false;
    if (brute_force) {
        // Naive hit algorithm
        float min_dist = 1e5f;
        for (const auto &triangle : triangles_) {
            HitRecord curr_hit_record;
            if (triangle.Hit(ray, &curr_hit_record)) {
                if (curr_hit_record.distance < min_dist) {
                    *hit_record = curr_hit_record;
                    min_dist = curr_hit_record.distance;
                }
            }
        }
        if (min_dist + 1.0 < 1e5f) {
            hit_record->material = material_;
            return true;
        }
        return false;
    } else {
        bool ret = OctreeHit(root_, ray, hit_record);
        if (ret) {
            hit_record->material = material_;
        }
        return ret;
    }
}

bool Mesh::IsFaceInsideBox(const std::vector<size_t>& face, const Point& bbox_min, const Point& bbox_max) const {
    for (size_t idx : face) {
        const auto& pt = vertices_[idx];
        for (int i = 0; i < 3; i++) {
            if (pt[i] < bbox_min[i] + 1e-6f) return false;
            if (pt[i] > bbox_max[i] - 1e-6f) return false;
        }
    }
    return true;
}

bool Mesh::IsRayIntersectBox(const Ray& ray, const Point& bbox_min, const Point& bbox_max) const {
    float t_min = -1e5f;
    float t_max =  1e5f;

    for (int i = 0; i < 3; i++) {
        if (glm::abs(ray.d[i]) < 1e-6f) {
            if (ray.o[i] < bbox_min[i] + 1e-6f || ray.o[i] > bbox_max[i] - 1e-6f) {
                t_min =  1e5f;
                t_max = -1e5f;
            }
        }
        else {
            if (ray.d[i] > 0.f) {
                t_min = glm::max(t_min, (bbox_min[i] - ray.o[i]) / ray.d[i]);
                t_max = glm::min(t_max, (bbox_max[i] - ray.o[i]) / ray.d[i]);
            }
            else {
                t_min = glm::max(t_min, (bbox_max[i] - ray.o[i]) / ray.d[i]);
                t_max = glm::min(t_max, (bbox_min[i] - ray.o[i]) / ray.d[i]);
            }
        }
    }

    return t_min + 1e-6f < t_max;
}

void Mesh::InsertFace(OctreeNode* u, size_t face_idx) {
    const Point& bbox_min = u->bbox_min;
    const Point& bbox_max = u->bbox_max;

    Vec bias = bbox_max - bbox_min;
    Vec half_bias = bias * 0.5f;

    bool inside_childs = false;

    for (size_t a = 0; a < 2; a++) {
        for (size_t b = 0; b < 2; b++) {
            for (size_t c = 0; c < 2; c++) {
                size_t child_idx = ((a << 2) | (b << 1) | c);
                Point curr_bbox_min = bbox_min + half_bias * Vec(float(a), float(b), float(c));
                Point curr_bbox_max = curr_bbox_min + half_bias;
                if (IsFaceInsideBox(f_ind_[face_idx], curr_bbox_min, curr_bbox_max)) {
                    if (u->childs[child_idx] == nullptr) {
                        tree_nodes_.emplace_back(new OctreeNode());
                        OctreeNode* child = tree_nodes_.back().get();
                        u->childs[child_idx] = tree_nodes_.back().get();
                        child->bbox_min = curr_bbox_min;
                        child->bbox_max = curr_bbox_max;
                    }
                    InsertFace(u->childs[child_idx], face_idx);
                    inside_childs = true;
                }
            }
        }
    }

    if (!inside_childs) {
        u->face_index.push_back(face_idx);
    }
}

bool Mesh::OctreeHit(OctreeNode* u, const Ray& ray, HitRecord* hit_record) const {
    if (!IsRayIntersectBox(ray, u->bbox_min, u->bbox_max)) {
        return false;
    }
    float distance = 1e5f;
    for (const auto& face_idx : u->face_index) {
        HitRecord curr_hit_record;
        if (triangles_[face_idx].Hit(ray, &curr_hit_record)) {
            if (curr_hit_record.distance < distance) {
                distance = curr_hit_record.distance;
                *hit_record = curr_hit_record;
            }
        }
    }

    for (const auto& child : u->childs) {
        if (child == nullptr) {
            continue;
        }
        HitRecord curr_hit_record;
        if (OctreeHit(child, ray, &curr_hit_record)) {
            if (curr_hit_record.distance < distance) {
                distance = curr_hit_record.distance;
                *hit_record = curr_hit_record;
            }
        }
    }
    return distance + 1 < 1e5f;
}


// Hittable list
void HittableList::PushHittable(const Hittable& hittable) {
    hittable_list_.push_back(&hittable);
}

bool HittableList::Hit(const Ray& ray, HitRecord *hit_record) const {
    float min_dist = 1e5f;
    for (const auto &hittable : hittable_list_) {
        HitRecord curr_hit_record;
        if (hittable->Hit(ray, &curr_hit_record)) {
            if (curr_hit_record.distance < min_dist) {
                *hit_record = curr_hit_record;
                min_dist = curr_hit_record.distance;
            }
        }
    }
    return min_dist + 1.0 < 1e4f;
}
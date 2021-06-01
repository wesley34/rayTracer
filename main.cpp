#include <iostream>
#include <vector>
#include "Common.h"
#include "Scene.h"
#include "Camera.h"
#include "Material.h"
#include "Hittable.h"
#include "Utils/lodepng.h"

const int kMaxTraceDepth = 5;

bool absoluteToleranceCompare(float x, float y)
{
    return std::fabs(x - y) <= std::numeric_limits<float>::epsilon() ;
}
Color TraceRay(const Ray& ray,
               const std::vector<LightSource>& light_sources,
               const Hittable& scene,
               int trace_depth);

Color Shade(const std::vector<LightSource>& light_sources,
            const Hittable& hittable_collection,
            const HitRecord& hit_record,
            int trace_depth) {
    // TODO: Add your code here.



    Color color;
    Color r_color;
    Color ambient_color;
    glm::vec3 shadow_ray;
    glm::vec3 reflected_ray = hit_record.reflection;
    glm::vec3 refracted_shadow_ray;
    glm::vec3 norm = hit_record.normal;
    glm::vec3 unit_norm;
    glm::vec3 norm_shadow_ray;
    glm::vec3 light_in_direction = hit_record.in_direction;
    glm::vec3 object_position =  hit_record.position;
    float reflection;
    bool is_hit = false;
    bool is_light = false;
    bool is_blocked = true;
    Material material = hit_record.material;
    HitRecord record_shadow ;
    glm::vec3 origin(0.0f,0.0f,0.0f);
 
 

    unit_norm = glm::normalize(norm);
   
    ambient_color = material.k_a * material.ambient;
    color = ambient_color;
   
    for (int i = 0; i < light_sources.size(); i++){
       
        shadow_ray = light_sources[i].position - object_position;
        float distance_between_obj_and_light = glm::length(shadow_ray);
        norm_shadow_ray = glm::normalize(shadow_ray);
      
        const Ray shadow_ray_casting(object_position+ 0.0001f*(unit_norm+norm_shadow_ray),norm_shadow_ray);
        
        is_hit = hittable_collection.Hit(shadow_ray_casting,&record_shadow);
        is_blocked = record_shadow.distance < distance_between_obj_and_light;


        if (glm::dot(norm_shadow_ray,unit_norm)>=0){
        
       
            if (!(is_hit && is_blocked)){
                 
                 // compute
            
                refracted_shadow_ray =  2 * glm::dot(norm_shadow_ray,unit_norm) * unit_norm -norm_shadow_ray ;
                reflection = glm::dot(refracted_shadow_ray,-light_in_direction);
                if (reflection > 0){
                    color += light_sources[i].intensity * (material.k_s * material.specular * abs(pow(reflection,material.sh)));  
                    color += light_sources[i].intensity * ( material.k_d * material.diffuse * glm::dot(norm_shadow_ray,unit_norm));  
                }
                else{
                    color += 1.0f * light_sources[i].intensity * ( material.k_d * material.diffuse * glm::dot(norm_shadow_ray,unit_norm));  
                }
               
            }

        }
        else{
            
            if (!(is_hit && is_blocked)){
                 
                 // compute
               
                refracted_shadow_ray =  2 * (glm::dot(norm_shadow_ray,unit_norm)+0.0f) * unit_norm -norm_shadow_ray ;
                reflection = glm::dot(refracted_shadow_ray,-light_in_direction);
                if (reflection > 0){
                    color += light_sources[i].intensity * (material.k_s * material.specular * abs(pow(reflection,material.sh)));  
                }
               
            }
            

        }
     
        
    }
    if (trace_depth < kMaxTraceDepth){
        if (material.k_s >= 0){
           if (glm::dot(unit_norm,glm::normalize(reflected_ray))>0){
            const Ray reflected (object_position+0.0001f*(unit_norm+glm::normalize(reflected_ray)),glm::normalize(reflected_ray));
            r_color= TraceRay(reflected, light_sources, hittable_collection, trace_depth++);
            color += material.k_s * r_color;
           }
      
           
           
        
        }
    }

    if (color[0] > 1){
        color[0] = 1;
    }
    if (color[1] > 1){
        color[1] = 1;
    }
    if (color[2] > 1){
        color[2] = 1;
    }

    


    return color;
}



Color TraceRay(const Ray& ray,
               const std::vector<LightSource>& light_sources,
               const Hittable& hittable_collection,
               int trace_depth) {
    // TODO: Add your code here.
    HitRecord record;
    Color color(0.0f, 0.0f, 0.0f);

    
    if (hittable_collection.Hit(ray,&record)){   
        return Shade(light_sources,hittable_collection,record,trace_depth++);
    }
  

    return color;
}



int main() {
    // TODO: Set your workdir (absolute path) here.
    const std::string work_dir("D://comp3271/Graphics_PA3_Release/Graphics_PA3_Release/");

    // Construct scene
    Scene scene(work_dir, "scene/spheres.toml");
    const Camera& camera = scene.camera_;
    int width = camera.width_;
    int height = camera.height_;

    std::vector<unsigned char> image(width * height * 4, 0);

    float progress = 0.f;

    // Traverse all pixels
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Color color(0.f, 0.f, 0.f);
            int count = 0;
            for (float bias_x = 0.25f; bias_x < 1.f; bias_x += .5f) {
                for (float bias_y = 0.25f; bias_y < 1.f; bias_y += .5f) {
                    Ray ray = camera.RayAt(float(x) + bias_x, float(y) + bias_y);
                    color += TraceRay(ray, scene.light_sources_, scene.hittable_collection_, 1);
                    count++;
                }
            }
            color /= float(count);
            int idx = 4 * ((height - y - 1) * width + x);
            for (int i = 0; i < 3; i++) {
                image[idx + i] = (uint8_t) (glm::min(color[i], 1.f - 1e-5f) * 256.f);
            }
            image[idx + 3] = 255;

            float curr_progress = float(x * height + y) / float(height * width);
            if (curr_progress > progress + 0.05f) {
                progress += 0.05f;
                std::cout << "Progress: " << progress << std::endl;
            }
        }
    }

    // Save result as png file
    std::vector<unsigned char> png;
    unsigned error = lodepng::encode(png, image, width, height);
    lodepng::save_file(png, work_dir + "output.png");
}


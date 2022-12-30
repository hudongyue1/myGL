# myGL
1. The simplest ray tracer is finished in v0.1, and it can be downloaded in https://github.com/hudongyue1/myGL/archive/refs/tags/lesson.zip
  
  The related post is [Scratchapixel: Simple Ray Tracing | HDY blog](https://hudongyue.com/2022/12/15/Scratchapixel-Simple-Ray-Tracing/)
  
  It can render the image like:

  ![v0.1](./resource/v0.1.png)

2. The ray tracer v0.2, and it can support to render the triangle mesh. It can be downloaded in https://github.com/hudongyue1/myGL/releases/tag/lesson1
  
  With a mesh model, it can triangulate it at first, then it can render it like:

  ![v0.2](./resource/v0.2.png)

Adding the choice for **Flat Shading** and **Smooth Shading**:

<div><table frame=void>	<!--用了<div>进行封装-->
	<tr>
        <td><div><center>	<!--每个格子内是图片加标题-->
        	<img src="./resource/Flat-Shading.png"
                 alt="Flat-Shading"
                 height="300"/>	<!--高度设置-->
        	Flat-Shading	<!--标题1-->
        </center></div></td>    
     	<td><div><center>	<!--第二张图片-->
    		<img src="./resource/Smooth-Shading.png"
                 alt="Smooth-Shading"
                 height="300"/>	
    		Smooth-Shading
        </center></div></td>
	</tr>
</table></div>
3. The ray tracer v0.3, add some effects for direct illumination. The related post is [Scratchapixel: Light and Shading | HDY blog](https://hudongyue.com/2022/12/23/Scratchapixel-Light-and-Shading/)

   * Distance light, point light and shadow

   <div><table frame=void>	<!--用了<div>进行封装-->
   	<tr>
           <td><div><center>	<!--每个格子内是图片加标题-->
           	<img src="./resource/shadow1.png"
                    alt="distanceShadow1"
                    height="300"/>	<!--高度设置-->
           	Shadow1	<!--标题1-->
           </center></div></td>    
        	<td><div><center>	<!--第二张图片-->
       		<img src="./resource/shadow2.png"
                    alt="distanceShadow2"
                    height="300"/>	
       		Shadow2
           </center></div></td>
   	</tr>
   </table></div>

   * Several light

     <img src="./resource/severalLight.png" alt="severalLight" style="zoom:50%;" />

   * Reflection

     <img src="./resource/reflection.png" alt="reflection" style="zoom:50%;" />

   * Reflection and refraction

     <div><table frame=void>	<!--用了<div>进行封装-->
     	<tr>
             <td><div><center>	<!--每个格子内是图片加标题-->
             	<img src="./resource/reflectionAndRefractionWithoutSmoothShading.png"
                      alt="distanceShadow1"
                      height="300"/>	<!--高度设置-->
             	Without Smooth Shading	<!--标题1-->
             </center></div></td>    
          	<td><div><center>	<!--第二张图片-->
         		<img src="./resource/reflectionAndRefractionWithSmoothShading.png"
                      alt="distanceShadow2"
                      height="300"/>	
         		With Smooth Shading
             </center></div></td>
     	</tr>
     </table></div>

     

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using WPoint = System.Windows.Point;
using System.Windows.Input;
using MathKernel.VectorAlgebra;
using ObjectModel.SceneGraph;

namespace ViewModelCore.OpticalSceneEntities.SceneView
{
    public class OpenGLRenderScene : OpenGLScene
    {
        protected RenderSceneObject renderTree;

        public void Render(RenderSceneObject rso)
        {
            RenderObject = rso;
            renderTree = rso;
        }

        public virtual SortedList<double, RenderSceneObject> HitScene(WPoint pos)
        {
            return HitScene(pos, renderTree);
        }

        protected SortedList<double, RenderSceneObject> HitScene(WPoint pos, RenderSceneObject root)
        {
            // create container for hitted objects
            var hittedRenderObjects = new SortedList<double, RenderSceneObject>();

            if (root == null) return hittedRenderObjects;

            // get world coordinates, create array of elements
            var near = new Point(GL.UnProject(pos.X, Height - pos.Y, 0));
            var far = new Point(GL.UnProject(pos.X, Height - pos.Y, 1));
            var dir = far - near;
            dir.Normalize();

            // process hit
            root.Hit(hittedRenderObjects, GL, (int)pos.X, (int)(Height - pos.Y), new Vector(near, dir));
            return hittedRenderObjects;
        }

    }
}

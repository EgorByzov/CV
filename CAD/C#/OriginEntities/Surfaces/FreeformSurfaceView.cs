using System;
using System.Collections.Generic;
using MathKernel.VectorAlgebra;
using ObjectModel.OpticalSceneEntities.OriginEntities.Surfaces;
using ViewModelCore.Common;

namespace ViewModelCore.OpticalSceneEntities.OriginEntities.Surfaces
{
    public class FreeformSurfaceView : SurfaceView
    {
        #region fields

        #endregion

        #region Constructors

        static FreeformSurfaceView()
        {
            InitRules();
        }

        public FreeformSurfaceView(FreeformSurface surface) : base(surface)
        {
            Surface = surface;
            InitFields();
        }

        #endregion

        #region Properties

        public FreeformSurface Surface { get; private set; }
        public int PointsCount { get; private set; }

        #endregion

        #region Private
        private static void InitRules()
        {
            
        }

        private void InitFields()
        {
            PointsCount = Surface.ControlPoints.Count;
        }
        #endregion

        #region Public

        #endregion
    }
}

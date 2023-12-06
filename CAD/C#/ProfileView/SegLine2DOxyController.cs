using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Common;
using MathKernel;
using MathKernel.VectorAlgebra;
using ObjectModel.OpticalSceneEntities.OriginEntities.Lines;
using OxyPlot;
using OxyPlot.Wpf;
using ViewModelCore.Controllers;

namespace ViewModelCore.OpticalSceneEntities.ProfileView
{
    internal sealed class SegLine2DOxyController : Line2DOxyController
    {        
        public SegLine2D SegLine { get; private set; }

        internal override Line2D Line2D { get { return SegLine; } }

        internal SegLine2DOxyController(SegLine2D line, CoordinateSystem castCS)
            : base(castCS)
        {
            SegLine = line;

            System.Windows.WeakEventManager<SegLine2D, ObjectChangedEventArgs>.AddHandler(SegLine, "ObjectChanged",
                (o, e) =>
                {
                    if (e.ID.Equals(Line2D.EventIDs.Geometry))
                        UpdateGeometry();
                });

            InitControlPoints();
            UpdateGeometry();
            UpdateLayout();
        }
    }
}

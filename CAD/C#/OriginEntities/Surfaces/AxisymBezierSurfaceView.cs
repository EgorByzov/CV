using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Input;
using ObjectModel.OpticalSceneEntities.OriginEntities.Surfaces;
using ViewModelCore.Common;

namespace ViewModelCore.OpticalSceneEntities.OriginEntities.Surfaces
{
    public class AxisymBezierSurfaceView : SurfaceView
    {
        private ICommand _updateSurfCommand;
        public AxisymBezierSurfaceView(AxisymBezierSurface surf)
            :base(surf)
        {
            SubscribeToModelObject<AxisymBezierSurface, AxisymBezierSurfaceView>(surf);
        }

        public ICommand UpdateSurfCommand
        {
            get
            {
                if (_updateSurfCommand == null)
                {
                    _updateSurfCommand = new Command(
                        param => UpdateSurface(),
                        param => true
                    );
                }
                return _updateSurfCommand;
            }
        }

        public void UpdateSurface()
        {
        }
    }
}

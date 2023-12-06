using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Controls.Primitives;
using Common;
using ViewModelCore.Controllers;

namespace ViewModelCore.OpticalSceneEntities.SceneView
{
    public class OpticalSceneTabControl : TabControl
    {
        #region Dependency
        public String GeneralHeader
        {
            get { return (String)GetValue(GeneralHeaderProperty); }
            set { SetValue(GeneralHeaderProperty, value); }
        }

        public static readonly DependencyProperty GeneralHeaderProperty =
            DependencyProperty.Register("GeneralHeader", typeof(String), typeof(OpticalSceneTabControl), new UIPropertyMetadata("General", OnHeaderChanged));

        public String SourcesHeader
        {
            get { return (String)GetValue(SourcesHeaderProperty); }
            set { SetValue(SourcesHeaderProperty, value); }
        }

        public static readonly DependencyProperty SourcesHeaderProperty =
            DependencyProperty.Register("SourcesHeader", typeof(String), typeof(OpticalSceneTabControl), new UIPropertyMetadata("Sources", OnHeaderChanged));

        public String RegistratorsHeader
        {
            get { return (String)GetValue(RegistratorsHeaderProperty); }
            set { SetValue(RegistratorsHeaderProperty, value); }
        }

        public static readonly DependencyProperty RegistratorsHeaderProperty =
            DependencyProperty.Register("RegistratorsHeader", typeof(String), typeof(OpticalSceneTabControl), new UIPropertyMetadata("Registrators", OnHeaderChanged));

        private static void OnHeaderChanged(DependencyObject d, DependencyPropertyChangedEventArgs e)
        {
            var uc = d as OpticalSceneTabControl;
            uc.UpdateHeaders();
        }

        #endregion Dependency

        private TabItem[] _items;
        private OSTreeView _osTreeView;
        private OSOneLevelTreeView _sourcesTreeView;
        private OSOneLevelTreeView _regsTreeView;
        private Project _curProject;

        public OpticalSceneTabControl()
        {
            _osTreeView = new OSTreeView(){BorderThickness = new Thickness(0)};
            _sourcesTreeView = new OSOneLevelTreeView() { BorderThickness = new Thickness(0) };
            _regsTreeView = new OSOneLevelTreeView() { BorderThickness = new Thickness(0) };

            _sourcesTreeView.PropertyChanged += (sender, args) => { if (args.PropertyName.Equals("NumElements")) UpdateHeaders(); };
            _regsTreeView.PropertyChanged += (sender, args) => { if (args.PropertyName.Equals("NumElements")) UpdateHeaders(); };

            _items = new TabItem[3] { new TabItem(), new TabItem(), new TabItem() };
            Items.Add(_items[0]); Items.Add(_items[1]); Items.Add(_items[2]);

            _items[0].Content = _osTreeView;
            _items[1].Content = _sourcesTreeView;
            _items[2].Content = _regsTreeView;

            InitCurrentProject();
            UpdateHeaders();


            ProjectController.StaticPropertyChanged += (sender, args) =>
            {
                if (args.PropertyName.Equals("CurrentProject"))
                    InitCurrentProject();
            };
            ApplicationPreferences.StaticPropertyChanged += (sender, args) =>
            {
                if (args.PropertyName.Equals("ColorPalette"))
                {
                    _osTreeView.UpdateColors();
                    _sourcesTreeView.UpdateColors();
                    _regsTreeView.UpdateColors();
                }
            };
        }

        private void InitCurrentProject()
        {
            _curProject = ProjectController.CurrentProject;

            _osTreeView.InitCurrentProject(OpticalSceneController.OpticalScenePresenter);
            _sourcesTreeView.InitCurrentProject(_curProject, ElementType.Source);
            _regsTreeView.InitCurrentProject(_curProject, ElementType.Registrator);
}

        private void UpdateHeaders()
        {
            if (_items.Length != 3) return;

            _items[0].Header = GeneralHeader;
            _items[1].Header = String.Format("{0} ({1})", SourcesHeader, _sourcesTreeView.Items.Count);
            _items[2].Header = String.Format("{0} ({1})", RegistratorsHeader, _regsTreeView.Items.Count);
        }
    }
}

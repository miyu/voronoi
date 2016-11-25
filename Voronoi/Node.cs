namespace Voronoi
{
    public class Node
    {
        private readonly Node[] _children = new Node[2];
        private Node _parent;

        public Node(object data)
        {
            Data = data;
        }

        public object Data { get; set; }

        public Node Parent { get { return _parent; } set { SetParent(value); } }

        private void SetParent(Node value)
        {
            if (_parent != value)
            {
                var oldParent = _parent;
                _parent = value;

                if (oldParent != null)
                {
                    if (oldParent.Left == this)
                    {
                        oldParent.Left = null;
                    }
                    else
                    {
                        oldParent.Right = null;
                    }
                }
            }
        }

        public Node Left
        {
            get { return GetChild(0); }
            set { SetChild(0, value); }
        }

        public Node Right
        {
            get { return GetChild(1); }
            set { SetChild(1, value); }
        }

        private Node GetChild(int index) => _children[index];
        private void SetChild(int index, Node node)
        {
            var existing = _children[index];
            if (existing?.Parent == this) existing.Parent = null;
            _children[index] = node;
            if (node != null) node.Parent = this;
        }
    }
}
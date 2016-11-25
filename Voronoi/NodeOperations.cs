using System;

namespace Voronoi
{
    public static class NodeOperations
    {
        public static bool IsLeaf(this Node root) => root.Left == null && root.Right == null;

        public static bool TryGetLeftAncestor(this Node root, out Node leftAncestor)
        {
            Node rightAncestor;
            FindDirectionalAncestors(root, out leftAncestor, out rightAncestor);
            return leftAncestor != null;
        }

        public static bool TryGetRightAncestor(this Node root, out Node rightAncestor)
        {
            Node leftAncestor;
            FindDirectionalAncestors(root, out leftAncestor, out rightAncestor);
            return rightAncestor != null;
        }

        public static void FindDirectionalAncestors(
            this Node root,
            out Node leftAncestor,
            out Node rightAncestor)
        {
            leftAncestor = rightAncestor = null;
            while (root.Parent != null && (leftAncestor == null || rightAncestor == null))
            {
                if (root.Parent.Left == root && rightAncestor == null)
                {
                    rightAncestor = root.Parent;
                } else if (root.Parent.Right == root && leftAncestor == null)
                {
                    leftAncestor = root.Parent;
                }
                root = root.Parent;
            }
        }

        public static void ValidateTree(this Node node)
        {
            if (node == null) return;

            if (node.Left != null)
            {
                if (node.Left.Parent != node) throw new Exception();
                ValidateTree(node.Left);
            }

            if (node.Right != null)
            {
                if (node.Right.Parent != node) throw new Exception();
                ValidateTree(node.Right);
            }

            if ((node.Left == null) != (node.Right == null)) throw new Exception();
        }

        public static Node GetLeftmost(this Node root)
        {
            while (root.Left != null) root = root.Left;
            return root;
        }

        public static Node GetRightmost(this Node root)
        {
            while (root.Right != null) root = root.Right;
            return root;
        }

        public static bool TryGetLeftLeaf(this Node leaf, out Node result)
        {
            if (!leaf.IsLeaf()) throw new ArgumentException();

            Node leftAncestor;
            if (leaf.TryGetLeftAncestor(out leftAncestor))
            {
                result = leftAncestor.Left.GetRightmost();
                return true;
            }
            result = null;
            return false;
        }

        public static bool TryGetRightLeaf(this Node leaf, out Node result)
        {
            if (!leaf.IsLeaf()) throw new ArgumentException();

            Node leftAncestor;
            if (leaf.TryGetRightAncestor(out leftAncestor))
            {
                result = leftAncestor.Right.GetLeftmost();
                return true;
            }
            result = null;
            return false;
        }

        public static void ReplaceNode(Node existing, Node replacement, ref Node root)
        {
            var existingParent = existing.Parent;
            if (existingParent == null)
            {
                root = replacement;
                replacement.Parent = null;
            }
            else if (existingParent.Left == existing)
            {
                existingParent.Left = replacement;
            }
            else
            {
                existingParent.Right = replacement;
            }
        }
    }
}
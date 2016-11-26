using System;
using System.Collections.Generic;

namespace Voronoi
{
    public class PriorityQueue<TItem>
    {
        private const int HeapBranchFactor = 2;
        private const int InitialCapacity = 1 + HeapBranchFactor;

        private readonly IComparer<TItem> _comparer;
        private TItem[] _storage = new TItem[InitialCapacity];
        private int _size;

        public PriorityQueue(IComparer<TItem> comparer)
        {
            _comparer = comparer;
        }

        public int Count => _size;

        public PriorityQueue<TItem> Copy()
        {
            var copy = new PriorityQueue<TItem>(_comparer);
            copy._storage = new TItem[_storage.Length];
            Array.Copy(_storage, copy._storage, _size);
            copy._size = _size;
            return copy;
        }

        public void Enqueue(TItem item)
        {
            // expand heap 1 level if storage is too small to fit another item.
            if (_size == _storage.Length)
            {
                var newStorage = new TItem[_storage.Length * HeapBranchFactor + 1];
                Array.Copy(_storage, newStorage, _storage.Length);
                _storage = newStorage;
            }

            // imaginary put of item into index _size
            _size++;

            // start percolating inserted item up from its current index, _size
            var childIndex = _size - 1;
            while (childIndex > 0)
            {
                // if inserted item is greater than or equal to its parent, stop percolation.
                var parentIndex = (childIndex - 1) / HeapBranchFactor;
                if (_comparer.Compare(item, _storage[parentIndex]) >= 0)
                    break;

                // Otherwise shift its parent down and continue percolating one level closer to root.
                _storage[childIndex] = _storage[parentIndex];
                childIndex = parentIndex;
            }

            // actually place item into the heap
            _storage[childIndex] = item;
        }

        public bool TryPeek(out TItem item)
        {
            item = _size > 0 ? _storage[0] : default(TItem);
            return _size > 0;
        }

        public bool TryDequeue(out TItem item) {
            if (_size == 0)
            {
                item = default(TItem);
                return false;
            }

            // remove heap head
            item = _storage[0];
            _size--;

            // remove heap tail
            var tail = _storage[_size];
            _storage[_size] = default(TItem);

            // percolate heap tail down from root; at iteration start, _storage[parentIndex] needs to be filled.
            var parentIndex = 0;
            while (true)
            {
                var childrenStartIndexInclusive = Math.Min(parentIndex * HeapBranchFactor + 1, _size);
                var childrenEndIndexExclusive = Math.Min(parentIndex * HeapBranchFactor + HeapBranchFactor + 1, _size);
                var childCount = childrenEndIndexExclusive - childrenStartIndexInclusive;

                // cannot percolate down beyond leaf node
                if (childCount == 0)
                {
                    _storage[parentIndex] = tail;
                    break;
                }

                // find smallest child
                var smallestChildIndex = childrenStartIndexInclusive;
                for (var i = childrenStartIndexInclusive + 1; i < childrenEndIndexExclusive; i++)
                {
                    if (_comparer.Compare(_storage[smallestChildIndex], _storage[i]) > 0)
                    {
                        smallestChildIndex = i;
                    }
                }

                // if smallest child is greater than the percolate down subject (heap tail), end percolation
                if (_comparer.Compare(_storage[smallestChildIndex], tail) > 0)
                {
                    _storage[parentIndex] = tail;
                    break;
                }

                // shift smallest child up the heap, continue percolation from its old index
                _storage[parentIndex] = _storage[smallestChildIndex];
                parentIndex = smallestChildIndex;
            }
            return true;
        }

        public TItem Dequeue()
        {
            TItem result;
            if (!TryDequeue(out result))
                throw new InvalidOperationException();
            return result;
        }

        public bool Any() => _size != 0;
    }
}
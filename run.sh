SCHEDULING_METHODS=(dynamic guided static)

CHUNK_SIZES=(1 10 100 1000)

THREAD_NUMS=(1 2 4 8)

M=(1000 10000 100000 10000000)

EXEC=$1

EXEC_OUTPUT_FILE=/dev/null

OUTPUT=""

EXEC_TIMES=()

for schedule in "${SCHEDULING_METHODS[@]}"; do

  for chunk in "${CHUNK_SIZES[@]}"; do
    export OMP_SCHEDULE="$schedule,$chunk"

    for m in "${M[@]}"; do
      OUTPUT="$m, $schedule, $chunk, "
      i=0
      for thread in "${THREAD_NUMS[@]}"; do
        export OMP_NUM_THREADS=$thread
        start=$(date +%s.%N)
        $EXEC "$m" >$EXEC_OUTPUT_FILE
        end=$(date +%s.%N)
        time=$(echo "$end - $start" | bc -l)
        EXEC_TIMES[$i]=$time
        OUTPUT="$OUTPUT, $time"
        i=$((i + 1))
      done

      for ((i = 1; i < $((${#THREAD_NUMS[@]})); i++)); do
        speed=$(echo "${EXEC_TIMES[0]} / ${EXEC_TIMES[$i]}" | bc -l)
        OUTPUT="$OUTPUT, $speed"
      done

      echo "$OUTPUT"
    done
  done
done

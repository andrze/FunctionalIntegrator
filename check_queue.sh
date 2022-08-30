export DISPLAY=:0
export DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus

catalog='/home/andrzej/Documents/Uczelnia/Anizotropie/FunctionalIntegrator/'
nq=$catalog.new_queue
oq=$catalog.old_queue

queue=$(ssh achlebicki@cluster.uy 'python ~/get_qstat.py') || return

echo "$queue" | head -n 3 > $nq

cmp --silent $oq $nq || notify-send "Zmiana w kolejce" "$(cat $nq)"

mv $nq $oq


